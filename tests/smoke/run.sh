#!/usr/bin/env bash
# poly-suite smoke test: bin/run.sh --reuse-prep (pgsc_calc scoring via Nextflow + Docker -> grade ->
# report) on a 19-record slice of a GIAB HG002 genotype-prepped VCF, with two small PGS Catalog scorefiles.
# Checks the raw score SUM and the matched-variant count per score against tests/smoke/expected.tsv and
# exits non-zero on any mismatch. PGS003442 equals the full-genome, panel-calibrated HG002 run exactly;
# PGS005161 differs from it by one site (chr8:142695299) that the calibrated run excluded at matching.
#
#   tests/smoke/run.sh [OUTDIR]        # default OUTDIR: scratch/smoke (gitignored)
#
# Needs: Docker, Nextflow (on PATH or $NF), network on first use (pulls pgsc_calc and its containers).
# Not covered: BAM input / genotype_prep.sh, the raw --vcf restriction path, ancestry calibration (--panel),
# the fastpath. Inputs: see tests/smoke/build_inputs.py. The checkout must sit under $OG_MOUNT (default /data).
set -euo pipefail
cd "$(dirname "$0")/../.."
S=tests/smoke
OUT=${1:-scratch/smoke}
export NF=${NF:-$(command -v nextflow)}
rm -rf "$OUT"; mkdir -p "$OUT/sfcache"
cp "$S"/scorefiles/*.txt.gz "$OUT/sfcache/"   # run.sh writes its target-list cache here; $OUT/scorefiles is its own staging dir
t0=$(date +%s)
bin/run.sh --sample HG002 --reuse-prep "$S/HG002.smoke.vcf.gz" --pgs PGS003442,PGS005161 --sex male \
  --scorefile-cache "$OUT/sfcache" --outdir "$OUT" > "$OUT/run.log" 2>&1 \
  || { echo "SMOKE FAIL: bin/run.sh exited non-zero (see $OUT/run.log)"; tail -5 "$OUT/run.log"; exit 1; }
echo "run.sh finished in $(( $(date +%s) - t0 )) s"

python3 - "$OUT" "$S/expected.tsv" <<'PY'
import csv, glob, gzip, sys
out, exp = sys.argv[1], sys.argv[2]
agg = glob.glob(f"{out}/score/**/aggregated_scores*.txt.gz", recursive=True)
sums = {}
with gzip.open(agg[0], "rt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        sums[r["PGS"].split("_")[0]] = float(r["SUM"])
graded = {r["pgs_id"]: r for r in csv.DictReader(open(f"{out}/score/pgs_scores.tsv"), delimiter="\t")}
bad = 0
for e in csv.DictReader(open(exp), delimiter="\t"):
    got_sum, got_n = sums.get(e["pgs_id"]), graded.get(e["pgs_id"], {}).get("n_matched")
    ok = got_sum is not None and abs(got_sum - float(e["sum"])) < 1e-6 and got_n == e["n_matched"]
    bad += not ok
    print(f"{'PASS' if ok else 'FAIL'} {e['pgs_id']}: SUM {got_sum} (expected {e['sum']}), "
          f"n_matched {got_n} (expected {e['n_matched']})")
sys.exit(1 if bad else 0)
PY
echo "SMOKE PASS"
