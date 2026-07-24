#!/usr/bin/env bash
# poly-suite end-to-end orchestrator — one entry point composing the stages:
#   [BAM]  harmonize scorefiles -> genotype-prep at scoring loci -> score -> grade -> report
#   [VCF]  score -> grade -> report
#   [reuse] --reuse-prep cached.vcf.gz  score -> grade -> report (skip harmonize+prep; prep-once)
# Optional ancestry calibration (--panel) and absolute risk (--sex, or inferred from --bam).
#
# usage:
#   run.sh --sample S --pgs "PGS...,PGS..." --outdir DIR \
#          (--vcf full_genotype.vcf.gz | --bam sample.bam --bootstrap-vcf boot.vcf.gz) \
#          [--panel panel.tar.zst] [--sex male|female] [--reuse-prep cached.vcf.gz]
#          [--work-cache DIR] [--ref REF] [--dry-run]
set -euo pipefail
cd "$(dirname "$0")/.."

SAMPLE= PGS= OUTDIR= VCF= BAM= PANEL= SEX= BOOT= DRY= REUSE_PREP= WORK_CACHE= SCOREFILE_CACHE= PGS_META=
BATCH= MAX_BATCH_BYTES=525000000  # --batch: size-aware bins cap peak RSS. 525MB -> ~35GB worst bin (measured: peak_rss=36.3*GB+15.8), safe 2-3 wide. See docs/memory-optimization-plan.md
REF=/data/alvin/ref/GRCh38/hg38.canonical.fa
NF=/home/alvin/bin/nextflow
while [ $# -gt 0 ]; do
  case "$1" in
    --sample) SAMPLE=$2; shift 2;; --pgs) PGS=$2; shift 2;;
    --outdir) OUTDIR=$2; shift 2;; --vcf) VCF=$2; shift 2;;
    --bam) BAM=$2; shift 2;; --panel) PANEL=$2; shift 2;;
    --sex) SEX=$2; shift 2;; --ref) REF=$2; shift 2;;
    --bootstrap-vcf) BOOT=$2; shift 2;; --reuse-prep) REUSE_PREP=$2; shift 2;;
    --work-cache) WORK_CACHE=$2; shift 2;; --scorefile-cache) SCOREFILE_CACHE=$2; shift 2;;
    --pgs-meta) PGS_META=$2; shift 2;;
    --batch) BATCH=1; shift;; --max-batch-bytes) MAX_BATCH_BYTES=$2; shift 2;;
    --dry-run) DRY=1; shift;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$SAMPLE" ] && [ -n "$PGS" ] && [ -n "$OUTDIR" ] || {
  echo "need --sample --pgs --outdir (and --vcf or --bam)"; exit 2; }
[ -n "$VCF" ] || [ -n "$BAM" ] || [ -n "$REUSE_PREP" ] || { echo "need --vcf, --bam, or --reuse-prep"; exit 2; }
mkdir -p "$OUTDIR"
# pgsc_calc resolves samplesheet path_prefix relative to its own base — all paths absolute
OUTDIR=$(readlink -f "$OUTDIR")
[ -n "$VCF" ] && VCF=$(readlink -f "$VCF")
[ -n "$BOOT" ] && BOOT=$(readlink -f "$BOOT")
[ -n "$PANEL" ] && PANEL=$(readlink -f "$PANEL")
[ -n "$REUSE_PREP" ] && REUSE_PREP=$(readlink -f "$REUSE_PREP")
[ -n "$PGS_META" ] && PGS_META=$(readlink -f "$PGS_META")

# --- scorefile cache is the DEFAULT path ------------------------------------------
# A persistent local cache of harmonized hmPOS scorefiles lets us skip pgsc_calc's
# DOWNLOAD_SCOREFILES step (the PGS-Catalog FTP fetch that hung ~66 min on 2026-07-23,
# and costs ~13 min even when healthy). If --scorefile-cache was not given, fall back
# to $POLY_SCOREFILE_CACHE, then to the first known persistent cache dir that exists.
# cache_scorefile_paths() is all-or-nothing: if the cache does not fully cover $PGS,
# the run transparently falls back to --pgs_id (download), so outputs never change.
if [ -z "$SCOREFILE_CACHE" ]; then
  for d in "${POLY_SCOREFILE_CACHE:-}" "results/launch70/scorefile_cache" "cache/scorefiles"; do
    [ -n "$d" ] && [ -d "$d" ] && { SCOREFILE_CACHE="$d"; break; }
  done
  [ -n "$SCOREFILE_CACHE" ] && echo "[cache] defaulting --scorefile-cache -> $SCOREFILE_CACHE" >&2
fi
[ -n "$SCOREFILE_CACHE" ] && SCOREFILE_CACHE=$(readlink -f "$SCOREFILE_CACHE")

ANC=(); [ -n "$PANEL" ] && ANC=(--run_ancestry "$PANEL")
FMT=vcf   # score-input format; reuse-prep upgrades to pfile when a .pgen exists (skips VCF->pgen)

# Resolve cached raw hmPOS scorefiles for a PGS list: prints one absolute path per id,
# all-or-nothing (returns nonzero and prints nothing on the first miss). A full hit lets us
# skip the harmonize nextflow run (feed scoring_targets.py directly) and the scoring download.
cache_scorefile_paths() {
  local dir=$1 ids=$2 id f out=()
  for id in ${ids//,/ }; do
    f="$dir/${id}_hmPOS_GRCh38.txt.gz"
    [ -f "$f" ] || return 1
    out+=("$(readlink -f "$f")")
  done
  printf '%s\n' "${out[@]}"
}

plan() { echo "  $*"; }
echo "== poly-suite run: $SAMPLE =="
plan "scores : $PGS"
plan "input  : ${BAM:+BAM $BAM}${VCF:+VCF $VCF}"
plan "ancestry: ${PANEL:-none (uncalibrated)}"

if [ -n "$REUSE_PREP" ]; then
  SCORE_IN="${REUSE_PREP%.vcf.gz}"; SCORE_IN="${SCORE_IN%.pgen}"
  [ -f "$SCORE_IN.pgen" ] && FMT=pfile
  plan "reuse cached prep: $SCORE_IN ($FMT)  (skip harmonize + genotype-prep${FMT:+; $FMT=no VCF conversion})"
  [ -n "$BAM" ] && [ -z "$SEX" ] && plan "infer sex from BAM"
elif [ -n "$BAM" ]; then
  SCORE_IN="$OUTDIR/${SAMPLE}.prepped"
  if [ -n "$SCOREFILE_CACHE" ] && cache_scorefile_paths "$SCOREFILE_CACHE" "$PGS" >/dev/null 2>&1; then
    plan "step 1: scorefile-cache hit -> skip harmonize (no --bootstrap-vcf needed)"
  else
    BOOT=${BOOT:-$VCF}
    [ -n "$BOOT" ] || { echo "BAM input needs --bootstrap-vcf (or --scorefile-cache with every id cached)"; exit 2; }
    plan "step 1: harmonize scorefiles (bootstrap $BOOT)"
  fi
  plan "step 2: genotype-prep BAM at scoring loci -> $OUTDIR/${SAMPLE}.prepped.vcf.gz"
  [ -z "$SEX" ] && plan "step 2b: infer sex from BAM"
else
  SCORE_IN="${VCF%.vcf.gz}"
fi
plan "step 3: pgsc_calc score${PANEL:+ + --run_ancestry}"
plan "step 4: grade -> contract + provenance + report.html  (sex=${SEX:-auto/none})"

if [ -n "$DRY" ]; then echo "== dry-run: plan only, nothing executed =="; exit 0; fi

# ---- execution ----
if [ -n "$REUSE_PREP" ]; then
  [ -f "$REUSE_PREP" ] || { echo "reuse-prep VCF not found: $REUSE_PREP"; exit 1; }
  [ -z "$SEX" ] && [ -n "$BAM" ] && SEX=$(python3 bin/infer_sex.py "$BAM" | sed -n 's/.*sex=\([a-z]*\).*/\1/p')
elif [ -n "$BAM" ]; then
  # Full scorefile-cache hit -> feed the cached hmPOS files straight to genotype-prep and
  # skip the harmonize nextflow run entirely (one fewer pgsc_calc cold-start on the BAM path).
  SF=()
  [ -n "$SCOREFILE_CACHE" ] && mapfile -t SF < <(cache_scorefile_paths "$SCOREFILE_CACHE" "$PGS")
  if [ "${#SF[@]}" -gt 0 ]; then
    echo "[harmonize] scorefile-cache hit (${#SF[@]} files) -> skipping harmonize nextflow run"
  else
    [ -n "$SCOREFILE_CACHE" ] && echo "[harmonize] scorefile-cache miss -> running harmonize" >&2
    printf 'sampleset,path_prefix,chrom,format\n%s,%s,,vcf\n' "$SAMPLE" "${BOOT%.vcf.gz}" > "$OUTDIR/hm.samplesheet.csv"
    rm -rf "$OUTDIR/hm" "$OUTDIR/work_hm"
    # The harmonize run is the only DOWNLOAD_SCOREFILES-triggering step left once the
    # scorefile cache is the default (the scoring run below uses --scorefile on a hit).
    # A stalled PGS-Catalog FTP connection makes DOWNLOAD_SCOREFILES hang with no exit
    # (observed ~66 min idle), and Nextflow's local executor does not enforce the
    # process `time` directive — so bound the whole harmonize run with GNU `timeout`
    # and a small retry budget. It fails LOUD instead of hanging indefinitely.
    HM_TIMEOUT=${HM_TIMEOUT:-1500}   # seconds per attempt (25 min); healthy fetch << this
    HM_RETRIES=${HM_RETRIES:-2}
    TO=""; command -v timeout >/dev/null 2>&1 && TO="timeout -k 30 ${HM_TIMEOUT}s"
    hm_rc=1
    for hm_try in $(seq 1 $((HM_RETRIES + 1))); do
      [ "$hm_try" -gt 1 ] && echo "[harmonize] attempt $hm_try/$((HM_RETRIES + 1)) after timeout/failure" >&2
      $TO env NXF_ANSI_LOG=false "$NF" run pgscatalog/pgsc_calc -profile docker \
        --input "$OUTDIR/hm.samplesheet.csv" --target_build GRCh38 --pgs_id "$PGS" --min_overlap 0.1 \
        -c conf/rootless.config --outdir "$OUTDIR/hm" -work-dir "$OUTDIR/work_hm" > "$OUTDIR/harmonize.log" 2>&1
      hm_rc=$?
      [ "$hm_rc" -eq 0 ] && break
      [ "$hm_rc" -eq 124 ] && echo "[harmonize] attempt $hm_try TIMED OUT after ${HM_TIMEOUT}s (likely DOWNLOAD_SCOREFILES hang)" >&2
    done
    [ "$hm_rc" -eq 0 ] || { echo "harmonize failed/timed out (rc=$hm_rc) — check $OUTDIR/harmonize.log"; exit 1; }
    mapfile -t SF < <(python3 -c "import glob,os;s={};[s.setdefault(os.path.basename(p),p) for p in glob.glob('$OUTDIR/work_hm/**/normalised_*_hmPOS_GRCh38.txt.gz',recursive=True)];print(chr(10).join(sorted(s.values())))")
    [ "${#SF[@]}" -ge 1 ] || { echo "no harmonized scorefiles — check $OUTDIR/harmonize.log"; exit 1; }
  fi
  NPROC=${NPROC:-24} bin/genotype_prep.sh "$BAM" "$REF" "$OUTDIR/${SAMPLE}.prepped.vcf.gz" "${SF[@]}"
  [ -z "$SEX" ] && SEX=$(python3 bin/infer_sex.py "$BAM" | sed -n 's/.*sex=\([a-z]*\).*/\1/p')
fi

printf 'sampleset,path_prefix,chrom,format\n%s,%s,,%s\n' "$SAMPLE" "$SCORE_IN" "$FMT" > "$OUTDIR/samplesheet.csv"
rm -rf "$OUTDIR/score"

# --scorefile-cache: use locally-cached harmonized scoring files (--scorefile) instead of
# --pgs_id, skipping the ~40-58s DOWNLOAD_SCOREFILES step. Falls back to --pgs_id on any miss.
SCORE_SEL=(--pgs_id "$PGS"); SFDIR=
if [ -n "$SCOREFILE_CACHE" ]; then
  CACHED_SF=(); mapfile -t CACHED_SF < <(cache_scorefile_paths "$SCOREFILE_CACHE" "$PGS")
  if [ "${#CACHED_SF[@]}" -gt 0 ]; then
    SFDIR="$OUTDIR/scorefiles"; rm -rf "$SFDIR"; mkdir -p "$SFDIR"
    for f in "${CACHED_SF[@]}"; do ln -sf "$f" "$SFDIR/"; done
    SCORE_SEL=(--scorefile "$SFDIR/*.txt.gz")
  else echo "scorefile-cache miss -> using --pgs_id" >&2; fi
fi
# --work-cache: persistent shared work dir + -resume, so the 16GB panel extraction
# (EXTRACT_DATABASE, keyed on the panel only) is reused across runs instead of re-unpacked.
WORKDIR="$OUTDIR/work"; RESUME=()
if [ -n "$WORK_CACHE" ]; then WORKDIR=$(readlink -f "$WORK_CACHE"); mkdir -p "$WORKDIR"; RESUME=(-resume)
else rm -rf "$OUTDIR/work"; fi
[ -n "$BATCH" ] && [ -z "$SFDIR" ] && echo "[batch] ignored: needs --scorefile-cache (nothing to partition on the --pgs_id path)" >&2
if [ -n "$BATCH" ] && [ -n "$SFDIR" ]; then
  # --batch: split the score set into size-aware bins so each pgsc_calc invocation's
  # peak RSS stays bounded (the full 2.4GB score set drove a 108GB peak -> global OOM).
  # Bins share ONE per-sample work dir with -resume so the score-independent ancestry
  # steps (panel extract + FRAPOSA projection) compute in bin 1 and cache-hit in 2..K.
  # merge_scores.py concatenates the disjoint per-bin score rows back into one score dir.
  mkdir -p "$WORKDIR"; rm -rf "$OUTDIR/score" "$OUTDIR/score_bins"
  mapfile -t BINS < <(python3 bin/partition_scores.py "$SFDIR" "$MAX_BATCH_BYTES" 2>>"$OUTDIR/score.log")
  echo "[batch] ${#BINS[@]} bins (budget $((MAX_BATCH_BYTES/1000000)) MB), shared work dir, -resume ancestry" >&2
  BIN_OUT=(); k=0
  for line in "${BINS[@]}"; do
    k=$((k+1)); bo="$OUTDIR/score_bins/bin_$k"; bsf="$OUTDIR/score_bins/sf_$k"
    mkdir -p "$bsf"; IFS=$'\t' read -ra FILES <<< "$line"
    for f in "${FILES[@]}"; do ln -sf "$f" "$bsf/"; done
    echo "[batch] bin $k/${#BINS[@]}: ${#FILES[@]} scores -> $bo" >&2
    NXF_ANSI_LOG=false "$NF" run pgscatalog/pgsc_calc -profile docker \
      --input "$OUTDIR/samplesheet.csv" --target_build GRCh38 --scorefile "$bsf/*.txt.gz" --min_overlap 0.1 "${ANC[@]}" \
      -c conf/rootless.config --outdir "$bo" -work-dir "$WORKDIR" -resume >> "$OUTDIR/score.log" 2>&1
    BIN_OUT+=("$bo")
  done
  mkdir -p "$OUTDIR/score"
  python3 bin/merge_scores.py "$OUTDIR/score" "$SAMPLE" "${BIN_OUT[@]}" >> "$OUTDIR/score.log" 2>&1
else
  NXF_ANSI_LOG=false "$NF" run pgscatalog/pgsc_calc -profile docker \
    --input "$OUTDIR/samplesheet.csv" --target_build GRCh38 "${SCORE_SEL[@]}" --min_overlap 0.1 "${ANC[@]}" \
    -c conf/rootless.config --outdir "$OUTDIR/score" -work-dir "$WORKDIR" "${RESUME[@]}" > "$OUTDIR/score.log" 2>&1
fi

[ -n "$SEX" ] && echo "$SEX" > "$OUTDIR/score/sample_sex.txt"
# grade_pgs reads pgs_catalog_meta.json from the score dir for trait names + evidence
# grades; drop the run's own select_pgs metadata there so custom-named metas resolve
# (without this, a non-default score set falls back to the starter meta -> trait==pgs_id).
[ -n "$PGS_META" ] && cp -f "$PGS_META" "$OUTDIR/score/pgs_catalog_meta.json"

# Neanderthal-ancestry %: a scoring-style dosage over an archaic-introgression SNP
# panel -> score/neanderthal.tsv (read by report_html's Deep-ancestry card). Runs
# BEFORE grade_pgs so its report render picks up the card. Best-effort / non-fatal.
#   NOTE: the shipped panel is a SEED (below the provisional threshold) -> reports a
#   PROVISIONAL % only. And the PGS-only prepped VCF does NOT cover the archaic loci,
#   so a real % needs the panel force-genotyped off the BAM. Production command:
#     bin/genotype_prep.sh <bam> <ref> <out.panel.vcf.gz>  # over a panel-derived scorefile
#     python3 bin/neanderthal.py "$OUTDIR/score" --targets <out.panel.vcf.gz>
NEAN_VCF=
[ -n "$VCF" ] && NEAN_VCF="$VCF"                                   # full-genotype VCF input
[ -z "$NEAN_VCF" ] && [ -n "$REUSE_PREP" ] && NEAN_VCF="$REUSE_PREP"
if [ -n "$NEAN_VCF" ] && [ -f "$NEAN_VCF" ]; then
  python3 bin/neanderthal.py "$OUTDIR/score" --vcf "$NEAN_VCF" --sample "$SAMPLE" \
    >>"$OUTDIR/score.log" 2>&1 || echo "note: neanderthal step skipped (see score.log)"
else
  python3 bin/neanderthal.py "$OUTDIR/score" --sample "$SAMPLE" \
    >>"$OUTDIR/score.log" 2>&1 || echo "note: neanderthal wiring stub skipped (see score.log)"
fi

python3 bin/grade_pgs.py "$OUTDIR/score"

# Reported discrimination (AUROC/C-index/R2) from the PGS Catalog -> run-local
# cache, then re-render so the report shows it. Bounded + best-effort: a slow or
# offline Catalog must never fail or hang a completed run (report already exists).
TO=""; command -v timeout >/dev/null 2>&1 && TO="timeout 180"
if $TO python3 bin/cache_performance.py "$OUTDIR/score/pgs_scores.tsv" \
        "$OUTDIR/score/pgs_performance.tsv" 2>>"$OUTDIR/score.log"; then
  python3 bin/report_html.py "$OUTDIR/score" >/dev/null
else
  echo "note: reported-performance fetch skipped (timeout/offline) — report omits 'reported perf'"
fi
echo "RUN_DONE -> $OUTDIR/score  (pgs_scores.tsv/.json + provenance.json + report.html)"
