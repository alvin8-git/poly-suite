#!/usr/bin/env python3
"""Rebuild the smoke-test inputs in tests/smoke/ from public GIAB data.

You do not need this to run the smoke test; the outputs are committed. It records how they were made.

  python3 tests/smoke/build_inputs.py PREPPED_VCF

PREPPED_VCF   GIAB HG002 all-sites VCF from `bin/genotype_prep.sh` (force-called at the PGS loci), run on
              the GIAB HG002 GRCh38 BAM.

Writes:
  HG002.smoke.vcf.gz(.tbi)  records at the scoring loci of the two scorefiles in scorefiles/; the
                            ##reference/##bcftools* header lines (local paths) are removed
"""
import glob, gzip, os, sys
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))

def loci():
    out = set()
    for sf in glob.glob(os.path.join(HERE, "scorefiles", "*_hmPOS_GRCh38.txt.gz")):
        with gzip.open(sf, "rt") as fh:
            hdr = None
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if hdr is None:
                    hdr = {k: i for i, k in enumerate(f)}
                    continue
                if f[hdr["hm_chr"]] and f[hdr["hm_pos"]]:
                    out.add(("chr" + f[hdr["hm_chr"]], int(f[hdr["hm_pos"]])))
    return sorted(out, key=lambda x: (int(x[0][3:]) if x[0][3:].isdigit() else 99, x[1]))

def main(vcf_in):
    sites = loci()
    # --- VCF slice -------------------------------------------------------------------------
    src = pysam.VariantFile(vcf_in)
    hdr = pysam.VariantHeader()
    for rec in src.header.records:
        if rec.key in ("reference",) or rec.key.startswith("bcftools"):
            continue
        hdr.add_record(rec)
    for s in src.header.samples:
        hdr.add_sample(s)
    out_vcf = os.path.join(HERE, "HG002.smoke.vcf.gz")
    dst = pysam.VariantFile(out_vcf, "wz", header=hdr)
    n = 0
    for c, p in sites:
        for rec in src.fetch(c, p - 1, p):
            if rec.pos == p:
                dst.write(rec); n += 1
    dst.close()
    pysam.tabix_index(out_vcf, preset="vcf", force=True)
    print(f"{n} VCF records at {len(sites)} loci -> {out_vcf}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(*sys.argv[1:])
