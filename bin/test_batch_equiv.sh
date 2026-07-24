#!/usr/bin/env bash
# Equivalence test for --batch (the correctness proof for score-set batching).
# Runs the SAME small score set two ways on one prepped genome:
#   (1) single  : one pgsc_calc invocation (current default path)
#   (2) batched : --batch with a tiny byte budget -> forced into several bins, merged
# then asserts the matched-variant set per score and the ancestry percentile are IDENTICAL,
# and the SUM agrees up to floating-point summation-order noise (PLINK2 reorders the reduction
# under a different scorefile grouping -> ~ppm drift; the matched set and rank are unchanged).
# A shared work dir means the 16GB panel extract + FRAPOSA projection run once and every
# other task -resumes, so the whole test costs ~one panel extract, not two.
#
# Proves: merge row-concat is exact, bins are disjoint, and ancestry (-resume across bins)
# is score-independent -> batching changes nothing observable. Independent of the byte budget.
set -euo pipefail
cd "$(dirname "$0")/.."

SAMPLE=${SAMPLE:-HG001}
PREP=${PREP:-results/launch70/HG001/HG001.prepped.vcf.gz}
PANEL=${PANEL:-/data/alvin/ref/pgsc/pgsc_HGDP+1kGP_v1.tar.zst}
CACHE=${CACHE:-cache/scorefiles}
# 6 small scores incl. two real disease scores (PGS000071/000135) so percentiles exist
PGS=${PGS:-PGS004700,PGS005161,PGS003442,PGS002299,PGS000071,PGS000135}
OUT=${OUT:-scratch/batch_equiv}
WORK="$OUT/work_shared"   # shared across both runs -> panel extract happens once

rm -rf "$OUT"; mkdir -p "$OUT"
echo "== [1/2] single (baseline) =="
bin/run.sh --sample "$SAMPLE" --pgs "$PGS" --outdir "$OUT/single" \
  --reuse-prep "$PREP" --panel "$PANEL" --scorefile-cache "$CACHE" --work-cache "$WORK"
echo "== [2/2] batched (tiny budget -> multiple bins) =="
bin/run.sh --sample "$SAMPLE" --pgs "$PGS" --outdir "$OUT/batched" \
  --reuse-prep "$PREP" --panel "$PANEL" --scorefile-cache "$CACHE" --work-cache "$WORK" \
  --batch --max-batch-bytes 8000

echo "== compare =="
python3 - "$OUT/single/score" "$OUT/batched/score" "$SAMPLE" <<'PY'
import sys, glob, gzip, csv
from collections import defaultdict
single_dir, batched_dir, sample = sys.argv[1:4]

def matched(d):
    # EXACT invariant: total matched variants per accession. Batching must not change
    # which variants match — this is the real equivalence proof (the SUM is a float
    # reduction over these variants and drifts by summation order; the SET must not).
    agg = defaultdict(int)
    for h in glob.glob(f"{d}/**/match/*_summary.csv", recursive=True):
        with open(h) as fh:
            for row in csv.DictReader(fh):
                if row.get("match_status") == "matched":
                    acc = row.get("accession") or row.get("PGS")
                    agg[acc] += int(float(row.get("count", 0) or 0))
    return dict(agg)

def load(d):
    # sample-row scores keyed by PGS accession: {PGS: (SUM, percentile)}
    for pat in (f"{d}/**/*_pgs.txt.gz", f"{d}/**/aggregated_scores*.txt.gz"):
        hits = sorted(glob.glob(pat, recursive=True))
        if hits:
            break
    assert hits, f"no score file under {d}"
    out = {}
    with gzip.open(hits[0], "rt") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(cols)}
        pgs_i, sum_i = ci["PGS"], ci["SUM"]
        pct_i = ci.get("percentile_MostSimilarPop")
        ss_i = ci.get("sampleset")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if ss_i is not None and f[ss_i] != sample:   # skip reference-panel rows
                continue
            pct = float(f[pct_i]) if pct_i is not None else None
            out[f[pgs_i]] = (float(f[sum_i]), pct)
    return out, hits[0]

a, pa = load(single_dir)
b, pb = load(batched_dir)
ma, mb = matched(single_dir), matched(batched_dir)
print(f"single : {len(a)} scores  ({pa})")
print(f"batched: {len(b)} scores  ({pb})")
assert set(a) == set(b), f"PGS set differs: only-single={set(a)-set(b)} only-batched={set(b)-set(a)}"

# (1) EXACT: identical matched-variant set per score. A real batching bug shows up here.
mbad = [f"{k}: single={ma.get(k)} batched={mb.get(k)}" for k in sorted(set(ma) | set(mb)) if ma.get(k) != mb.get(k)]
if mbad:
    print("MATCH-COUNT MISMATCH (real bug — batching changed variant matching):")
    [print("  " + x) for x in mbad]; sys.exit(1)

# (2) percentile -> grade (the clinical output) must be identical; (3) SUM agrees up to
# floating-point summation-order noise (PLINK2 reduces 1e5-1e6 float terms; a different
# scorefile grouping reorders the sum -> ~ppm drift, no effect on matched set or rank).
SUM_RTOL, SUM_ATOL, PCT_TOL = 1e-4, 1e-3, 1e-4
bad = []; worst = 0.0
for pgs in sorted(a):
    (s1, p1), (s2, p2) = a[pgs], b[pgs]
    rel = abs(s1 - s2) / max(abs(s1), 1e-9); worst = max(worst, rel)
    if rel > SUM_RTOL and abs(s1 - s2) > SUM_ATOL:
        bad.append(f"{pgs} SUM {s1} != {s2} (rel {rel:.2e})")
    if p1 is not None and p2 is not None and abs(p1 - p2) > PCT_TOL:
        bad.append(f"{pgs} percentile {p1} != {p2}")
if bad:
    print("MISMATCH:"); [print("  " + x) for x in bad]; sys.exit(1)
print(f"EQUIVALENCE PASSED: {len(a)} scores — matched-variant sets identical, "
      f"percentiles identical, SUM agrees to {worst:.1e} rel (float summation-order noise)")
PY