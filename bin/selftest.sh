#!/usr/bin/env bash
# poly-suite self-test: all module self-checks + unit tests. No network/BAM/panel.
# Run: bin/selftest.sh
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== module self-checks =="
python3 bin/absolute_risk.py
python3 bin/consensus.py

echo "== unit tests =="
python3 bin/tests.py

echo "== grader end-to-end on synthetic calibrated fixture =="
T=$(mktemp -d)/score; mkdir -p "$T"
python3 - "$T" <<'PY'
import gzip, sys
rows = ["sampleset\tIID\tPGS\tSUM\tDENOM\tAVG\tpercentile_MostSimilarPop\tMostSimilarPop",
        "S\tS\tPGS000004\t1.2\t500\t0.0024\t87.0\tEUR"]
with gzip.open(f"{sys.argv[1]}/aggregated_scores.txt.gz", "wt") as f:
    f.write("\n".join(rows) + "\n")
PY
python3 bin/grade_pgs.py "$(dirname "$T")" >/dev/null
for out in pgs_scores.tsv pgs_scores.json provenance.json report.html; do
  [ -s "$(dirname "$T")/$out" ] || { echo "FAIL: $out not produced"; exit 1; }
done
echo "  grader emits contract + provenance + html OK"
python3 bin/validate_contract.py "$(dirname "$T")" >/dev/null
echo "  contract schema + gate invariants OK"

echo "ALL SELFTESTS PASSED"
