#!/usr/bin/env python3
"""Infer sample sex from a BAM via samtools idxstats — normalized X/Y coverage
relative to the autosomes. Needed for sex-dimorphic absolute risk (e.g. CAD,
whose baseline incidence differs by sex). Fast (reads the index only).

Female: chrX ~1.0x autosome, chrY ~0.   Male: chrX ~0.5x, chrY ~0.5x.

Importable (infer_sex) and runnable:
  python3 bin/infer_sex.py BAM [OUT_DIR]   # writes OUT_DIR/sample_sex.txt
"""
import sys, os, subprocess


def infer_sex(bam):
    out = subprocess.run(["samtools", "idxstats", bam],
                         capture_output=True, text=True, timeout=180).stdout
    cov = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) >= 3 and int(f[1]) > 0:
            cov[f[0]] = int(f[2]) / int(f[1])           # mapped reads / contig length

    def c(name):
        return cov.get(name) or cov.get(name.replace("chr", "")) or 0.0

    autos = [c(f"chr{i}") for i in range(1, 23)]
    autos = [a for a in autos if a > 0]
    if not autos:
        return "unknown", {}
    auto = sum(autos) / len(autos)
    xr, yr = c("chrX") / auto, c("chrY") / auto
    if xr >= 0.80 and yr <= 0.15:
        sex = "female"
    elif xr <= 0.70 and yr >= 0.12:
        sex = "male"
    else:
        sex = "unknown"
    return sex, {"x_ratio": round(xr, 3), "y_ratio": round(yr, 3),
                 "autosome_cov": round(auto, 4)}


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: infer_sex.py BAM [OUT_DIR]")
    sex, m = infer_sex(sys.argv[1])
    print(f"sex={sex}  x_ratio={m.get('x_ratio')}  y_ratio={m.get('y_ratio')}  "
          f"autosome_cov={m.get('autosome_cov')}")
    if len(sys.argv) > 2:
        os.makedirs(sys.argv[2], exist_ok=True)
        p = os.path.join(sys.argv[2], "sample_sex.txt")
        with open(p, "w") as f:
            f.write(sex + "\n")
        print("wrote", p)
