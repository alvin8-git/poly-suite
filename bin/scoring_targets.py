#!/usr/bin/env python3
"""Build a bcftools `-C alleles` targets file (CHROM POS REF,ALT[,ALT2...]) from
PGS harmonized scorefiles, for force-genotyping a BAM at scoring loci.

v2: handles SNPs, indels (anchored VCF-style, e.g. C/CAAA, CT/C), and multiallelic
sites (multiple ALTs per position). REF/ALT are oriented against the reference
FASTA — the reference-matching allele is REF; for indels the longest ref-matching
allele wins (so a deletion's longer allele is correctly REF). Adds the 'chr'
prefix (scorefiles are Ensembl-style '1'; the BAM/ref are 'chr1'). Sites whose
alleles don't match the reference are dropped and reported.

usage: scoring_targets.py OUT.tsv REF.fa harmonized_scorefile.txt.gz...
"""
import sys, gzip, pysam

OUT, REF = sys.argv[1], sys.argv[2]
SCOREFILES = sys.argv[3:]
fa = pysam.FastaFile(REF)
contigs = set(fa.references)
order = {c: i for i, c in enumerate(fa.references)}
ACGT = set("ACGT")
# Autosomes only: chrX needs sex-aware (haploid-male) dosage that pgsc_calc/plink2
# reject without sample sex; chrX is a small coverage fraction. Documented limitation.
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def orient(chrom, pos, a1, a2):
    """(ref, alt): ref = the reference-matching allele (longest match, so an
    indel's deleted-sequence allele is REF); None if neither matches the genome."""
    L = max(len(a1), len(a2))
    ref_seq = fa.fetch(chrom, pos - 1, pos - 1 + L).upper()
    matches = [a for a in (a1, a2) if ref_seq[:len(a)] == a]
    if not matches:
        return None
    ref = max(matches, key=len)
    alt = a2 if ref == a1 else a1
    return None if ref == alt else (ref, alt)


sites = {}  # (chrom,pos) -> {"ref": ref, "alts": [alt,...]}
orient_cache = {}  # (chrom,pos,ea,oa) -> orient() result; overlapping scorefiles share loci,
                   # so this collapses N*variants FASTA fetches down to ~unique loci.
n_indel = skip_contig = skip_nonacgt = skip_refmismatch = skip_conflict = skip_sex = 0
for sf in SCOREFILES:
    with gzip.open(sf, "rt") as fh:
        # raw PGS-Catalog *_hmPOS_*.txt.gz carry a '#' comment header; normalised_*
        # files don't. Skip comments, then take the first real line as the header.
        for header in fh:
            if not header.startswith("#"):
                break
        else:
            continue                                    # empty / all-comment file
        ix = {c: i for i, c in enumerate(header.rstrip("\n").split("\t"))}
        # raw hmPOS files hold GRCh38 coords in hm_chr/hm_pos (chr_name/chr_position are
        # the ORIGINAL build); normalised_* files already put GRCh38 in chr_name/position.
        chrom_col = "hm_chr" if "hm_chr" in ix else "chr_name"
        pos_col = "hm_pos" if "hm_pos" in ix else "chr_position"
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                chrom = f[ix[chrom_col]]
            except IndexError:
                continue
            chrom = chrom if chrom.startswith("chr") else "chr" + chrom
            if chrom not in contigs:
                skip_contig += 1
                continue
            if chrom not in AUTOSOMES:      # exclude chrX/Y/M (sex-aware handling deferred)
                skip_sex += 1
                continue
            try:
                pos = int(f[ix[pos_col]])
            except (ValueError, KeyError, IndexError):
                continue                                # unmapped variant (blank hm_pos)
            ea, oa = f[ix["effect_allele"]].upper(), f[ix["other_allele"]].upper()
            if not oa and "hm_inferOtherAllele" in ix:  # raw hmPOS may leave other_allele blank
                oa = f[ix["hm_inferOtherAllele"]].upper()
            if not ea or not oa or not (set(ea) <= ACGT and set(oa) <= ACGT):
                skip_nonacgt += 1
                continue
            ck = (chrom, pos, ea, oa)
            if ck in orient_cache:
                o = orient_cache[ck]
            else:
                o = orient_cache[ck] = orient(chrom, pos, ea, oa)
            if o is None:
                skip_refmismatch += 1
                continue
            ref, alt = o
            if len(ref) != 1 or len(alt) != 1:
                n_indel += 1
            key = (chrom, pos)
            if key not in sites:
                sites[key] = {"ref": ref, "alts": [alt]}
            elif sites[key]["ref"] == ref:
                if alt not in sites[key]["alts"]:
                    sites[key]["alts"].append(alt)          # multiallelic
            else:
                skip_conflict += 1                          # conflicting REF orientation

rows = sorted(sites.items(), key=lambda kv: (order.get(kv[0][0], 1 << 30), kv[0][1]))
with open(OUT, "w") as out:
    for (chrom, pos), v in rows:
        out.write(f"{chrom}\t{pos}\t{v['ref']},{','.join(v['alts'])}\n")

sys.stderr.write(
    f"[scoring_targets] {len(rows):,} sites ({n_indel:,} indels) | dropped: "
    f"refmismatch={skip_refmismatch:,} off-contig={skip_contig:,} "
    f"non-ACGT={skip_nonacgt:,} ref-conflict={skip_conflict:,} sex-chr={skip_sex:,}\n")
