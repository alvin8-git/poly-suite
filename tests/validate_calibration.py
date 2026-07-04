#!/usr/bin/env python3
"""Validation harness (spec §5) — calibration self-consistency from the reference
panel that pgsc_calc --run_ancestry emits alongside the sample.

By construction, ancestry-normalized percentiles for the ~3,300 reference-panel
samples must be ~uniform on 0-100 (each sample sits at its own rank within its
most-similar population). A non-uniform reference distribution means the FRAPOSA
per-ancestry normalization is off for that score. This checks it per PGS — the
defensibility pillar the space review flagged as missing from raw PGS tools.

Importable (uniformity) + CLI: validate_calibration.py <_pgs.txt.gz | results_dir>
"""
import sys, os, gzip, glob, csv
from collections import defaultdict

THRESH = 0.06   # a decile off by >6 percentage points from 10% => suspicious


def find_pgs_file(arg):
    if os.path.isfile(arg):
        return arg
    hits = sorted(glob.glob(f"{arg}/**/*_pgs.txt.gz", recursive=True))
    return hits[0] if hits else None


def uniformity(percentiles):
    """(max |decile_fraction - 0.1|, deciles) over 0-100; 0 = perfectly uniform.
    None if too few points to judge."""
    n = len(percentiles)
    if n < 100:
        return None
    dec = [0] * 10
    for p in percentiles:
        dec[min(9, max(0, int(p // 10)))] += 1
    fracs = [c / n for c in dec]
    return max(abs(f - 0.1) for f in fracs), fracs


def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else "results"
    f = find_pgs_file(arg)
    if not f:
        sys.exit(f"no calibrated *_pgs.txt.gz under {arg} (run pgsc_calc --run_ancestry)")
    ref = defaultdict(list)
    with gzip.open(f, "rt") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        pcol = next((c for c in r.fieldnames if "percentile" in c.lower()), None)
        for row in r:
            if (row.get("sampleset") or "").lower() == "reference":
                try:
                    ref[row["PGS"]].append(float(row[pcol]))
                except (ValueError, KeyError, TypeError):
                    pass
    print(f"[validate_calibration] {os.path.relpath(f)}  "
          f"({len(ref)} scores, {sum(len(v) for v in ref.values()):,} reference rows)\n")
    bad = 0
    for pgs in sorted(ref):
        res = uniformity(ref[pgs])
        if res is None:
            print(f"  {pgs[:22]:24} n={len(ref[pgs])} (too few to test)")
            continue
        maxdev, _ = res
        mean = sum(ref[pgs]) / len(ref[pgs])
        ok = maxdev <= THRESH and 45 <= mean <= 55
        bad += 0 if ok else 1
        print(f"  {pgs[:22]:24} n={len(ref[pgs]):>5} mean={mean:5.1f} "
              f"max_decile_dev={maxdev:.3f}  {'OK' if ok else 'FLAG'}")
    print(f"\n{'ALL CALIBRATIONS SELF-CONSISTENT' if bad == 0 else f'{bad} score(s) FLAGGED'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
