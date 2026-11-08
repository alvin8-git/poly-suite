#!/usr/bin/env python3
"""Raw-VCF input helpers for run.sh --vcf.

pgsc_calc's PLINK2_VCF step passes no sex to plink2, so any VCF containing chrX
aborts ("chrX is present ... no sex information was provided"). run.sh therefore
imports autosomes only and uses the sex inferred here for absolute risk.

  python3 bin/vcf_input.py sex   < GT-lines   # one GT per line (bcftools query -f '[%GT]\n')
  python3 bin/vcf_input.py restrict-expr      # prints the bcftools -t exclusion expression
  python3 bin/vcf_input.py targets SCOREFILE_CACHE REF ID,ID,...  # prints 'path<TAB>n'; cached CHROM/POS .tsv.gz+.tbi
  python3 bin/vcf_input.py fill ALLELES.tsv.gz CHROM COUNTS_OUT  # stdin VCF text -> stdout: fill ALT of hom-ref target rows
  python3 bin/vcf_input.py afrac TARGETS.tsv.gz   # stdin: bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n'; explicit-GT fraction

Target restriction: pgsc_calc's PLINK2_VCF runs out of memory on a whole-genome all-sites VCF
(HG001: 2.88 G rows), so run.sh restricts every --vcf input to the PGS target positions
(bin/scoring_targets.py, the same list genotype_prep.sh builds) with `bcftools view -R`
before import. A record kept at a target position keeps its genotype, so hom-ref 0/0 rows of
an all-sites input survive; a variant-only input yields a sparse file (see explicit_fraction).

ALT fill: `bcftools call -m` (no -C alleles) writes hom-ref rows with ALT ".", and pgsc_calc's MATCH
needs REF and ALT to equal the scorefile's alleles, so every hom-ref target site would be "unmatched"
(THAL1 launch76: median match_rate 0.479 vs ~0.90 on the BAM path, where `call -C alleles` puts the
target ALT on hom-ref rows). fill_homref_alt() therefore sets ALT on hom-ref rows (ALT ".") whose
(CHROM, POS, REF) is in the allele table (<targets>.alleles.tsv.gz = CHROM POS REF ALT from
scoring_targets.py, the same builder genotype_prep.sh uses): GT stays 0/0 and AD/ADF/ADR gain a
trailing 0. Genuine variant rows are untouched. Limitation: one ALT per site - at a multi-ALT target
site the FIRST table ALT is used. It streams per contig, so it runs inside the parallel restriction.

infer_sex_from_gts(): chrX non-PAR heterozygosity, het/(het+homalt):
  >= 0.25 female, < 0.10 male, else unknown (caller falls back to autosomes-only
  and reports no sex for absolute risk).
"""
import gzip
import hashlib
import os
import subprocess
import sys

# GRCh38 chrX PAR1 = 10,001-2,781,479 ; PAR2 = 155,701,383-156,030,895
NONPAR_X = "chrX:2781480-155701382"
EXCLUDE = "^chrX,chrY,chrM,X,Y,M,MT"
MIN_SITES = 50  # fewer informative sites than this -> unknown


def classify_gt(gt):
    """'het' | 'homalt' | None for one GT string (diploid, phased or not)."""
    a = gt.replace("|", "/").split("/")
    if len(a) != 2 or "." in a:
        return None
    if a[0] != a[1]:
        return "het"
    return "homalt" if a[0] != "0" else None


def infer_sex_from_gts(gts, min_sites=MIN_SITES):
    het = hom = 0
    for g in gts:
        c = classify_gt(g.strip())
        if c == "het":
            het += 1
        elif c == "homalt":
            hom += 1
    n = het + hom
    if n < min_sites:
        return "unknown", {"het": het, "homalt": hom, "het_frac": None}
    f = het / n
    sex = "female" if f >= 0.25 else "male" if f < 0.10 else "unknown"
    return sex, {"het": het, "homalt": hom, "het_frac": round(f, 4)}


HERE = os.path.dirname(os.path.abspath(__file__))
MIN_EXPLICIT = 0.5  # below this the input is not an all-sites VCF


def targets_key(ids, ref):
    """Cache key: sorted PGS ids + reference + the target-builder's own source."""
    h = hashlib.sha1()
    h.update(",".join(sorted(ids)).encode())
    h.update(b"\0" + os.path.basename(ref).encode())
    with open(os.path.join(HERE, "scoring_targets.py"), "rb") as fh:
        h.update(b"\0" + fh.read())
    return h.hexdigest()[:12]


def alleles_path(targets_path):
    """CHROM POS REF ALT table (bgzip+tabix) that sits next to a CHROM/POS targets file."""
    return targets_path[:-len(".tsv.gz")] + ".alleles.tsv.gz"


def get_targets(cache_dir, ids, scorefile_dir, ref, log=sys.stderr):
    """Return (path, n_positions) of the bgzipped, tabix-indexed CHROM/POS list for `ids`.

    Built with scoring_targets.py (as genotype_prep.sh does) from
    <scorefile_dir>/<id>_hmPOS_GRCh38.txt.gz and cached as
    <cache_dir>/vcf_targets_<hash>.tsv.gz (+ .tbi, .n); a repeat call is a cache hit."""
    import pysam
    ids = sorted(ids)
    base = os.path.join(cache_dir, f"vcf_targets_{targets_key(ids, ref)}")
    out, cnt = base + ".tsv.gz", base + ".n"
    aout = alleles_path(out)
    if all(os.path.exists(x) for x in (out, out + ".tbi", cnt, aout, aout + ".tbi")):
        n = int(open(cnt).read())
        log.write(f"[targets] cache hit {out} ({n:,} positions)\n")
        return out, n
    sfs = [os.path.join(scorefile_dir, f"{i}_hmPOS_GRCh38.txt.gz") for i in ids]
    missing = [f for f in sfs if not os.path.isfile(f)]
    if missing:
        raise FileNotFoundError(f"scorefile cache lacks {len(missing)} id(s), e.g. {missing[0]}")
    tmp = f"{base}.tmp{os.getpid()}"
    full = tmp + ".full.tsv"
    subprocess.run([sys.executable, os.path.join(HERE, "scoring_targets.py"), full, ref] + sfs,
                   check=True, stdout=sys.stderr)
    n = 0
    with open(full) as fi, open(tmp + ".tsv", "w") as fo, open(tmp + ".alleles.tsv", "w") as fa:
        for line in fi:
            c, p, al = line.rstrip("\n").split("\t", 2)
            fo.write(f"{c}\t{p}\n")
            ref, alts = al.split(",", 1)
            fa.write(f"{c}\t{p}\t{ref}\t{alts.split(',')[0]}\n")   # first ALT per site
            n += 1
    os.remove(full)
    pysam.tabix_compress(tmp + ".alleles.tsv", tmp + ".alleles.tsv.gz", force=True)
    os.remove(tmp + ".alleles.tsv")
    pysam.tabix_index(tmp + ".alleles.tsv.gz", seq_col=0, start_col=1, end_col=1, zerobased=False,
                      force=True, keep_original=True)
    os.replace(tmp + ".alleles.tsv.gz.tbi", aout + ".tbi")
    os.replace(tmp + ".alleles.tsv.gz", aout)
    pysam.tabix_compress(tmp + ".tsv", tmp + ".tsv.gz", force=True)
    os.remove(tmp + ".tsv")
    pysam.tabix_index(tmp + ".tsv.gz", seq_col=0, start_col=1, end_col=1, zerobased=False,
                      force=True, keep_original=True)
    open(cnt, "w").write(f"{n}\n")
    os.replace(tmp + ".tsv.gz.tbi", out + ".tbi")
    os.replace(tmp + ".tsv.gz", out)
    log.write(f"[targets] built {out} ({n:,} positions)\n")
    return out, n


R_FORMAT_KEYS = ("AD", "ADF", "ADR")   # Number=R: one more (zero) value for the added ALT


def load_alleles(alleles_gz, chrom):
    """{pos: (ref, alt)} for one contig from the tabix-indexed allele table."""
    import pysam
    tb = pysam.TabixFile(alleles_gz)
    try:
        if chrom not in tb.contigs:
            return {}
        return {int(f[1]): (f[2], f[3]) for f in (l.split("\t") for l in tb.fetch(chrom))}
    finally:
        tb.close()


def fill_homref_alt(lines, alleles):
    """Generator over VCF text lines: set ALT on hom-ref rows (ALT '.') from `alleles`
    ({pos: (ref, alt)} for the row's contig). Returns lines; counts in fill_homref_alt.stats
    (homref rows seen, rows filled). A row is filled only when its REF equals the table REF."""
    st = fill_homref_alt.stats = {"homref": 0, "filled": 0}
    for ln in lines:
        if ln.startswith("#"):
            yield ln
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 10 or f[4] != ".":
            yield ln
            continue
        st["homref"] += 1
        hit = alleles.get(int(f[1]))
        if not hit or hit[0] != f[3].upper():
            yield ln
            continue
        f[4] = hit[1]
        keys = f[8].split(":")
        idx = [i for i, k in enumerate(keys) if k in R_FORMAT_KEYS]
        if idx:
            for j in range(9, len(f)):
                v = f[j].split(":")
                for i in idx:
                    if i < len(v) and v[i] not in (".", ""):
                        v[i] += ",0"
                f[j] = ":".join(v)
        st["filled"] += 1
        yield "\t".join(f) + "\n"


def explicit_fraction(lines, target_rows):
    """lines: '%CHROM\\t%POS[\\t%GT]' rows of the restricted VCF; target_rows: iterable of
    (chrom, pos) target positions (strings/ints).

    The denominator is the targets inside the input's own footprint (per contig, first..last
    record position of the restricted file), so a region slice is judged on the region it
    covers and a whole-genome input on the whole genome. Returns
    (n_explicit, n_targets_in_footprint, fraction): hom-ref 0/0 counts as explicit, ./. does not."""
    seen, span = set(), {}
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        pos = int(f[1])
        lo, hi = span.get(f[0], (pos, pos))
        span[f[0]] = (min(lo, pos), max(hi, pos))
        if len(f) > 2 and any(classify_any(g) for g in f[2:]):
            seen.add((f[0], pos))
    n = sum(1 for c, p in target_rows if c in span and span[c][0] <= int(p) <= span[c][1])
    return len(seen), n, (len(seen) / n if n else 0.0)


def classify_any(gt):
    """True when GT is a called genotype (no missing allele), hom-ref included."""
    a = gt.replace("|", "/").split("/")
    return bool(a) and "." not in a and gt != ""


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "sex":
        sex, m = infer_sex_from_gts(sys.stdin)
        print(f"sex={sex}  chrX_nonPAR het={m['het']} homalt={m['homalt']} het_frac={m['het_frac']}")
    elif cmd == "nonpar-region":
        print(NONPAR_X)
    elif cmd == "restrict-expr":
        print(EXCLUDE)
    elif cmd == "targets":
        path, n = get_targets(sys.argv[2], sys.argv[4].split(","), sys.argv[2], sys.argv[3])
        print(f"{path}\t{n}")
    elif cmd == "fill":
        # fill ALLELES.tsv.gz CHROM COUNTS_OUT ; stdin VCF text -> stdout
        al = load_alleles(sys.argv[2], sys.argv[3])
        for out in fill_homref_alt(sys.stdin, al):
            sys.stdout.write(out)
        st = fill_homref_alt.stats
        with open(sys.argv[4], "w") as fh:
            fh.write(f"{st['homref']}\t{st['filled']}\n")
    elif cmd == "afrac":
        with gzip.open(sys.argv[2], "rt") as tf:
            k, n, fr = explicit_fraction(sys.stdin, (l.split("\t")[:2] for l in tf))
        print(f"explicit={k} targets_in_footprint={n} fraction={fr:.4f}")
    else:
        sys.exit(__doc__)
