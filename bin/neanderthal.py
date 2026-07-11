#!/usr/bin/env python3
"""poly-suite Neanderthal-ancestry estimation (standalone, scoring-style).

A dosage over an archaic-introgression-informative SNP panel: at each panel locus
the sample's count of the Neanderthal-derived (archaic) allele is summed, divided
by alleles genotyped, and mapped to a genome-wide Neanderthal-ancestry percentage
via a reference calibration. This is NOT a PGS Catalog score — the panel and the
calibration live in resources/ and the loci are force-genotyped off the BAM the
same way the PGS scoring loci are (bin/genotype_prep.sh), because the PGS-only
prepped VCF does not cover the archaic tag SNPs.

Contract (what OmniGen consumes):
  <score_dir>/neanderthal.tsv   columns: sample, neanderthal_pct, method

The panel that ships (resources/neanderthal_panel.tsv) is a documented SEED — the
schema plus a handful of well-established introgression loci — so percentages from
it are PROVISIONAL and flagged as such in the `method` column until the full
curated tag-SNP set is dropped in (recipe in the panel header). The wiring and the
dosage/calibration math are real and unit-tested; only the panel size is seed-scale.

Usage:
  neanderthal.py <score_dir> --vcf full_genotype.vcf.gz [--sample S]
  neanderthal.py <score_dir> --targets <targets.vcf.gz>   (force-genotyped panel VCF)
"""
import os, sys, csv, json, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PANEL = os.path.join(HERE, "..", "resources", "neanderthal_panel.tsv")
DEFAULT_CALIB = os.path.join(HERE, "..", "resources", "neanderthal_calibration.json")

CONTRACT_COLS = ["sample", "neanderthal_pct", "method"]
MIN_PANEL_SITES = 500          # below this the % is provisional (seed panel), not for interpretation
MIN_GENOTYPED = 50             # below this we cannot estimate at all


# --------------------------------------------------------------------------- panel / calibration
def load_panel(path=DEFAULT_PANEL):
    """[(chrom, pos, ref, archaic_allele, weight), ...] from a TSV with a '#'-comment
    header. Columns: chrom pos ref archaic_allele [weight]. 'chr' prefix normalized on."""
    out = []
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = [c.strip().lower() for c in f]
                continue
            row = dict(zip(header, f))
            chrom = row["chrom"]
            chrom = chrom if chrom.startswith("chr") else "chr" + chrom
            try:
                pos = int(row["pos"])
            except (KeyError, ValueError):
                continue
            ref = (row.get("ref") or "").upper()
            arch = (row.get("archaic_allele") or row.get("archaic") or "").upper()
            if not ref or not arch:
                continue
            try:
                w = float(row.get("weight", 1) or 1)
            except ValueError:
                w = 1.0
            out.append((chrom, pos, ref, arch, w))
    return out


def load_calibration(path=DEFAULT_CALIB):
    try:
        with open(path) as fh:
            return json.load(fh)
    except OSError:
        return {"slope": 1.0, "intercept": 0.0, "clamp_max": 10.0,
                "note": "default (calibration file missing)"}


# --------------------------------------------------------------------------- core math (pure, tested)
def archaic_fraction(archaic_count, total_alleles):
    """Fraction of panel alleles that are the archaic (Neanderthal-derived) allele."""
    if not total_alleles:
        return None
    return archaic_count / total_alleles


def pct_from_counts(archaic_count, total_alleles, calib=None):
    """Map archaic-allele fraction -> genome-wide Neanderthal-ancestry %.
    Linear calibration (pct = intercept + slope*fraction), clamped to [0, clamp_max].
    Returns None when nothing was genotyped."""
    frac = archaic_fraction(archaic_count, total_alleles)
    if frac is None:
        return None
    c = calib or load_calibration()
    pct = c.get("intercept", 0.0) + c.get("slope", 1.0) * frac
    lo, hi = 0.0, c.get("clamp_max", 10.0)
    return round(max(lo, min(hi, pct)), 2)


def method_string(n_panel, n_genotyped, calib=None):
    c = calib or {}
    tag = c.get("version", "v1")
    prov = " PROVISIONAL-SEED (expand panel before interpreting)" if n_panel < MIN_PANEL_SITES else ""
    return (f"dosage-over-archaic-tag-panel-{tag} "
            f"({n_genotyped}/{n_panel} sites genotyped, calibrated to 1kGP){prov}")


# --------------------------------------------------------------------------- genotype extraction
def archaic_counts_from_vcf(vcf_path, panel, sample=None):
    """(archaic_count, total_alleles, n_sites_genotyped) by reading the sample's
    genotypes at panel loci. Needs pysam + a VCF that covers the panel (a full-genotype
    VCF, or a panel-targeted force-genotyped VCF from genotype_prep.sh)."""
    import pysam
    vf = pysam.VariantFile(vcf_path)
    samp = sample or (list(vf.header.samples)[0] if vf.header.samples else None)
    by_locus = {}
    for (chrom, pos, ref, arch, w) in panel:
        by_locus.setdefault((chrom, pos), (ref, arch, w))
        by_locus.setdefault((chrom.replace("chr", ""), pos), (ref, arch, w))
    ac = tot = ngt = 0
    for (chrom, pos), (ref, arch, w) in by_locus.items():
        try:
            recs = list(vf.fetch(chrom, pos - 1, pos))
        except (ValueError, OSError):
            continue
        for rec in recs:
            if rec.pos != pos:
                continue
            call = rec.samples.get(samp) if samp else None
            if call is None or None in (call.get("GT") or (None,)):
                continue
            alleles = [rec.ref] + list(rec.alts or [])
            gt = call["GT"]
            counted = False
            for a in gt:
                if a is None or a >= len(alleles):
                    continue
                tot += 1
                counted = True
                if (alleles[a] or "").upper() == arch:
                    ac += 1
            if counted:
                ngt += 1
            break
    return ac, tot, ngt


# --------------------------------------------------------------------------- contract IO
def write_contract(score_dir, sample, pct, method, out=None):
    """Write <score_dir>/neanderthal.tsv (sample, neanderthal_pct, method)."""
    out = out or os.path.join(score_dir, "neanderthal.tsv")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CONTRACT_COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerow({"sample": sample,
                    "neanderthal_pct": ("NA" if pct is None else pct),
                    "method": method})
    return out


def read_contract(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --------------------------------------------------------------------------- CLI
def main(argv=None):
    ap = argparse.ArgumentParser(description="Neanderthal-ancestry % over an archaic SNP panel")
    ap.add_argument("score_dir", help="output score dir (writes neanderthal.tsv here)")
    ap.add_argument("--vcf", help="full-genotype VCF covering the panel loci")
    ap.add_argument("--targets", help="panel force-genotyped VCF (from genotype_prep.sh)")
    ap.add_argument("--panel", default=DEFAULT_PANEL)
    ap.add_argument("--calib", default=DEFAULT_CALIB)
    ap.add_argument("--sample", default=None)
    args = ap.parse_args(argv)

    panel = load_panel(args.panel)
    calib = load_calibration(args.calib)
    vcf = args.vcf or args.targets
    sample = args.sample or os.path.basename(os.path.dirname(os.path.abspath(args.score_dir))) or "sample"

    if not vcf:
        method = method_string(len(panel), 0, calib) + " | NO VCF SUPPLIED — wiring only"
        write_contract(args.score_dir, sample, None, method)
        print(f"[neanderthal] no --vcf/--targets given; wrote stub contract "
              f"(panel={len(panel)} sites). Real % needs a force-genotyped panel VCF.")
        return 0

    ac, tot, ngt = archaic_counts_from_vcf(vcf, panel, args.sample)
    if ngt < MIN_GENOTYPED:
        method = method_string(len(panel), ngt, calib) + " | INSUFFICIENT genotyped sites"
        write_contract(args.score_dir, sample, None, method)
        print(f"[neanderthal] only {ngt} panel sites genotyped (<{MIN_GENOTYPED}) — reporting NA")
        return 0
    pct = pct_from_counts(ac, tot, calib)
    method = method_string(len(panel), ngt, calib)
    out = write_contract(args.score_dir, sample, pct, method)
    print(f"[neanderthal] {sample}: {pct}% (archaic {ac}/{tot} alleles over {ngt} sites) -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
