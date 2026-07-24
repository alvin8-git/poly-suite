#!/usr/bin/env bash
# Measurement pass (memory-optimization plan §4): fit the scorefile-bytes -> peak_rss
# curve so --max-batch-bytes B is measured, not a proxy. Runs pgsc_calc scoring on
# HG001's prepped VCF over increasing scorefile subsets, WITHOUT --run_ancestry (the
# FORMAT/MATCH memory hogs are score-driven, ancestry is not — skipping it makes each
# run cheap and isolates the driver we care about).
#
# SERIAL by construction (one subset at a time) and it WAITS for any in-flight nextflow
# to finish first, so it never collides with the equivalence test -> no OOM. Each subset
# gets a FRESH work dir with NO -resume, so FORMAT/MATCH actually re-run and produce a
# fresh peak_rss in the trace; the work dir is pruned after the trace is copied out.
set -euo pipefail
cd "$(dirname "$0")/.."

NF=${NF:-/home/alvin/bin/nextflow}
SAMPLE=${SAMPLE:-HG001}
PREP=${PREP:-results/launch70/HG001/HG001.prepped.vcf.gz}
CACHE=${CACHE:-cache/scorefiles}
OUT=${OUT:-results/mem_measure}
WORKROOT=${WORKROOT:-work/mem_measure}

SCORE_IN="${PREP%.vcf.gz}"; FMT=vcf
[ -f "$SCORE_IN.pgen" ] && { SCORE_IN="${SCORE_IN%.pgen}"; FMT=pfile; }
PREP_ABS=$(readlink -f "$PREP"); SCORE_IN=$(readlink -f "$SCORE_IN")
[ -f "$PREP_ABS" ] || { echo "prepped input not found: $PREP"; exit 1; }

mkdir -p "$OUT" "$WORKROOT"
: > "$OUT/subsets.tsv"
printf 'sampleset,path_prefix,chrom,format\n%s,%s,,%s\n' "$SAMPLE" "$SCORE_IN" "$FMT" > "$OUT/samplesheet.csv"

# Don't collide with the equivalence test (or any other run). Wait it out.
while pgrep -f 'test_batch_equiv.sh' >/dev/null || pgrep -x nextflow >/dev/null; do
  echo "[wait] in-flight nextflow/equiv-test running; holding $(date '+%H:%M:%S')"
  sleep 120
done
echo "[go] no in-flight nextflow — starting measurement $(date '+%H:%M:%S')"

# Rank scorefiles by size desc; build increasing subsets.
mapfile -t RANKED < <(ls -1 "$CACHE"/*.txt.gz | xargs -I{} stat -c '%s %n' {} | sort -rn | awk '{print $2}')
big8=("${RANKED[@]:0:8}")
small=("${RANKED[@]:8}")   # the 128-file tail

# subset name -> file list (increasing total bytes; small tail last as the high anchor)
run_subset() {
  local name=$1; shift
  local files=("$@")
  local sfdir="$OUT/sf_$name"; rm -rf "$sfdir"; mkdir -p "$sfdir"
  local bytes=0
  for f in "${files[@]}"; do ln -sf "$(readlink -f "$f")" "$sfdir/"; bytes=$((bytes + $(stat -c '%s' "$f"))); done
  printf '%s\t%d\t%d\n' "$name" "$bytes" "${#files[@]}" >> "$OUT/subsets.tsv"

  local wd="$WORKROOT/$name"; rm -rf "$wd"; mkdir -p "$wd"
  local od="$OUT/out_$name"; rm -rf "$od"
  echo "[run] $name: ${#files[@]} files, $((bytes/1000000)) MB  $(date '+%H:%M:%S')"
  # NO --run_ancestry, NO -resume, fresh work dir -> forces FORMAT/MATCH to execute + be traced.
  NXF_ANSI_LOG=false "$NF" run pgscatalog/pgsc_calc -profile docker \
    --input "$OUT/samplesheet.csv" --target_build GRCh38 --scorefile "$sfdir/*.txt.gz" --min_overlap 0.1 \
    -c conf/rootless.config --outdir "$od" -work-dir "$wd" > "$OUT/$name.log" 2>&1 \
    || { echo "[fail] $name — see $OUT/$name.log (keeping work dir for triage)"; return 0; }

  local tr; tr=$(find "$od" -name 'execution_trace*.txt' 2>/dev/null | sort | tail -1)
  if [ -n "$tr" ]; then cp "$tr" "$OUT/$name.trace.txt"; echo "[ok] $name -> $OUT/$name.trace.txt"; fi
  rm -rf "$wd"   # work dir is pure cache once the trace is out
}

run_subset b1    "${big8[@]:0:1}"
run_subset b2    "${big8[@]:0:2}"
run_subset b4    "${big8[@]:0:4}"
run_subset b8    "${big8[@]}"
run_subset small "${small[@]}"

echo "== fit =="
python3 bin/fit_budget.py "$OUT" | tee "$OUT/budget.txt"
echo "[done] $(date '+%H:%M:%S')  -> $OUT/budget.txt"
