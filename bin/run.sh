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
#          [--check-input]   # --vcf only: restrict to PGS targets + plink2 import test, then exit
#          [--nf-timeout SEC]  # wall limit per pgsc_calc scoring run (default NF_TIMEOUT=28800 = 8 h)
#          [--fastpath [--fastpath-cache DIR]]  # fast per-sample scoring (bin/fastpath.py; env POLY_FASTPATH=1).
#               Same scores as pgsc_calc, sample-independent work cached once (`bin/fastpath.py build`).
#               Needs --panel and a complete cache; a missing/mismatching cache is a hard error (no fallback).
#               Env: POLY_FASTPATH_THREADS (default NPROC or 8), POLY_FASTPATH_MEM_GB (default 100; server needs more).
#               pgsc_calc stays the default until the fastpath is accepted (docs/fastpath.md).
# --vcf: pgsc_calc's PLINK2_VCF runs out of memory on a whole-genome all-sites VCF and passes no
#   sex to plink2 (chrX aborts the run). The VCF is first restricted to the PGS target positions
#   (parallel per-contig `bcftools view -r chrN -T targets`) -> $OUTDIR/<S>.targeted.vcf.gz, then staged to
#   $OUTDIR/<S>.autosomes.vcf.gz (chrX/chrY/chrM excluded); sex (--sex, else inferred from
#   chrX non-PAR het) is kept for absolute risk. Input must be an ALL-SITES VCF. See README.
set -euo pipefail
cd "$(dirname "$0")/.."

SAMPLE= PGS= OUTDIR= VCF= BAM= PANEL= SEX= BOOT= DRY= REUSE_PREP= WORK_CACHE= SCOREFILE_CACHE= PGS_META=
FASTPATH=${POLY_FASTPATH:-} FASTPATH_CACHE=
CHECK= NF_TIMEOUT=${NF_TIMEOUT:-28800} NF_FAIL_GRACE=${NF_FAIL_GRACE:-600}
BATCH= MAX_BATCH_BYTES=525000000  # --batch: size-aware bins cap peak RSS. 525MB -> ~35GB worst bin (measured: peak_rss=36.3*GB+15.8), safe 2-3 wide. See docs/memory-optimization-plan.md
REF=/data/alvin/ref/GRCh38/hg38.canonical.fa
NF=${NF:-/home/alvin/bin/nextflow}
while [ $# -gt 0 ]; do
  case "$1" in
    --sample) SAMPLE=$2; shift 2;; --pgs) PGS=$2; shift 2;;
    --outdir) OUTDIR=$2; shift 2;; --vcf) VCF=$2; shift 2;;
    --bam) BAM=$2; shift 2;; --panel) PANEL=$2; shift 2;;
    --sex) SEX=$2; shift 2;; --ref) REF=$2; shift 2;;
    --bootstrap-vcf) BOOT=$2; shift 2;; --reuse-prep) REUSE_PREP=$2; shift 2;;
    --work-cache) WORK_CACHE=$2; shift 2;; --scorefile-cache) SCOREFILE_CACHE=$2; shift 2;;
    --pgs-meta) PGS_META=$2; shift 2;;
    --fastpath) FASTPATH=1; shift;; --fastpath-cache) FASTPATH_CACHE=$2; FASTPATH=1; shift 2;;
    --batch) BATCH=1; shift;; --max-batch-bytes) MAX_BATCH_BYTES=$2; shift 2;;
    --dry-run) DRY=1; shift;; --check-input) CHECK=1; shift;;
    --nf-timeout) NF_TIMEOUT=$2; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$CHECK" ] && [ -z "$PGS" ] && { echo "--check-input needs --pgs (the target positions come from the PGS set)"; exit 2; }
[ -n "$CHECK" ] && [ -z "$VCF" ] && { echo "--check-input needs --vcf"; exit 2; }
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
  SCORE_IN="$OUTDIR/${SAMPLE}.autosomes"
  plan "step 1a: index input if needed; restrict to PGS target positions (per-contig bcftools view -r -T, parallel) -> $OUTDIR/${SAMPLE}.targeted.vcf.gz"
  plan "step 1: stage targeted VCF -> $SCORE_IN.vcf.gz (chrX/chrY/chrM excluded: pgsc_calc passes plink2 no sex)"
  [ -z "$SEX" ] && plan "step 1b: infer sex from VCF chrX non-PAR het; unknown -> no sex"
fi
if [ -n "$FASTPATH" ]; then plan "step 3: FASTPATH score (bin/fastpath.py; cached match/format/intersect)"; else plan "step 3: pgsc_calc score${PANEL:+ + --run_ancestry}"; fi
plan "step 4: grade -> contract + provenance + report.html  (sex=${SEX:-auto/none})"

if [ -n "$DRY" ]; then echo "== dry-run: plan only, nothing executed =="; exit 0; fi

# ---- execution ----
BCF_IMG=${BCFTOOLS_IMG:-quay.io/biocontainers/bcftools:1.21--h8b25389_0}
PLINK2_IMG=${PLINK2_IMG:-ghcr.io/pgscatalog/plink2:2.00a5.10}
bcf() { docker run --rm -v /data:/data -w "$PWD" "$BCF_IMG" bcftools "$@"; }

# (check-input mirrors pgsc_calc's PLINK2_VCF flags: --new-id-max-allele-len 100 missing, --chr 1-22,X,Y,XY)
# Run one pgsc_calc/nextflow invocation: own process group, wall-clock limit ($NF_TIMEOUT),
# and a watchdog -- once the log reports a failed pipeline ("Execution cancelled" / "Pipeline
# completed with errors") and nextflow has not exited within $NF_FAIL_GRACE s, kill the group
# (Nextflow otherwise waits for unrelated running tasks for hours). Returns non-zero on failure.
# usage: run_nf LOGFILE nextflow-args...   (log is appended)
run_nf() {
  local log=$1; shift
  local start=$SECONDS failed_at= pid pgid rc=0
  : >> "$log"
  NXF_ANSI_LOG=false setsid "$NF" "$@" >> "$log" 2>&1 &
  pid=$!; pgid=$pid
  while kill -0 "$pid" 2>/dev/null; do
    sleep 10
    if [ $((SECONDS - start)) -ge "$NF_TIMEOUT" ]; then
      echo "[nextflow] TIMEOUT after ${NF_TIMEOUT}s -> killing process group $pgid ($log)" >&2
      kill -TERM -- -"$pgid" 2>/dev/null; sleep 30; kill -KILL -- -"$pgid" 2>/dev/null; wait "$pid" 2>/dev/null
      return 124
    fi
    if [ -z "$failed_at" ] && grep -qE "Execution cancelled|Pipeline completed with errors|^ERROR ~" "$log" 2>/dev/null; then
      failed_at=$SECONDS; echo "[nextflow] pipeline reported an error; allowing ${NF_FAIL_GRACE}s to exit ($log)" >&2
    fi
    if [ -n "$failed_at" ] && [ $((SECONDS - failed_at)) -ge "$NF_FAIL_GRACE" ]; then
      echo "[nextflow] still running ${NF_FAIL_GRACE}s after error -> killing process group $pgid" >&2
      kill -TERM -- -"$pgid" 2>/dev/null; sleep 30; kill -KILL -- -"$pgid" 2>/dev/null; wait "$pid" 2>/dev/null
      return 1
    fi
  done
  wait "$pid" || rc=$?
  return $rc
}
if [ -n "$REUSE_PREP" ]; then
  [ -f "$REUSE_PREP" ] || { echo "reuse-prep VCF not found: $REUSE_PREP"; exit 1; }
  [ -z "$SEX" ] && [ -n "$BAM" ] && SEX=$(python3 bin/infer_sex.py "$BAM" | sed -n 's/.*sex=\([a-z]*\).*/\1/p')
elif [ -z "$BAM" ] && [ -n "$VCF" ]; then
  # Stage the raw VCF for pgsc_calc: autosomes only (see header). Sex for absolute risk is
  # --sex, else inferred from the VCF's chrX non-PAR heterozygosity (indexed VCFs only).
  [ -f "$VCF" ] || { echo "VCF not found: $VCF"; exit 1; }
  THREADS=${THREADS:-8}
  T0=$SECONDS
  # (1) index is mandatory: -R and the chrX sex inference both seek through it.
  if [ -f "$VCF.tbi" ] || [ -f "$VCF.csi" ]; then IDX=present
  else
    t=$SECONDS; echo "[vcf] no .tbi/.csi next to $VCF: indexing (bcftools index -t)"
    bcf index -t --threads "$THREADS" "$VCF" || { echo "cannot index $VCF (is its directory writable?)"; exit 1; }
    echo "[vcf] indexed in $((SECONDS - t)) s"; IDX="created ($((SECONDS - t)) s)"
  fi
  IN_ROWS=$(bcf index -n "$VCF" 2>/dev/null || echo "?")
  case "$(bcf index -s "$VCF" 2>/dev/null | head -1)" in chr*|"") ;; *)
    echo "VCF contigs are not chr-prefixed (GRCh38 'chr1' style required for the PGS target list)"; exit 1;; esac
  # (2) sex from the full input (indexed chrX non-PAR query), before restriction.
  if [ -z "$SEX" ]; then
    SEXLINE=$(bcf query -r "$(python3 bin/vcf_input.py nonpar-region)" -i 'GT="alt"' -f '[%GT\n]' "$VCF" | python3 bin/vcf_input.py sex)
    echo "[vcf] $SEXLINE"
    SEX=$(echo "$SEXLINE" | sed -n 's/^sex=\([a-z]*\).*/\1/p')
    [ "$SEX" = unknown ] && SEX=
  fi
  # (3) PGS target positions (same list genotype_prep.sh builds), cached per PGS-set hash.
  [ -n "$SCOREFILE_CACHE" ] && cache_scorefile_paths "$SCOREFILE_CACHE" "$PGS" >/dev/null 2>&1 || {
    echo "--vcf needs a scorefile cache covering every --pgs id (to build the target positions); see --scorefile-cache"; exit 1; }
  TL=$(python3 bin/vcf_input.py targets "$SCOREFILE_CACHE" "$REF" "$PGS") || { echo "target list build failed"; exit 1; }
  TARGETS=${TL%%$'\t'*}; NTARGETS=${TL##*$'\t'}
  # (4) restrict through the index (-R), then drop sex chromosomes on the small file.
  t=$SECONDS
  # A single `-R` over ~10 M point regions behaves like a scan (>67 min on the 20 GB HG001 file),
  # so restrict per contig in parallel: `-r <chrom> -T <targets>` streams the targets inside one
  # region; parts are concatenated in contig order (`concat --naive`, identical headers).
  PAR=${NPROC:-24}   # parallel per-contig workers (same NPROC convention as genotype_prep.sh)
  PARTS="$OUTDIR/${SAMPLE}.targeted_parts"; rm -rf "$PARTS"; mkdir -p "$PARTS"
  mapfile -t CHROMS < <(bcf index -s "$VCF" | cut -f1 | grep -xE 'chr([1-9]|1[0-9]|2[0-2])' | sort -V)
  ALLELES=${TARGETS%.tsv.gz}.alleles.tsv.gz
  export BCF_IMG PWD VCF TARGETS PARTS ALLELES
  printf '%s\n' "${CHROMS[@]}" | xargs -P "$PAR" -I{} bash -c \
    'docker run --rm -v /data:/data -w "$PWD" "$BCF_IMG" bcftools view -r "$1" -T "$TARGETS" -Ov "$VCF" \
       | python3 bin/vcf_input.py fill "$ALLELES" "$1" "$PARTS/$1.fill" | bgzip -c > "$PARTS/$1.vcf.gz"; \
     [ "${PIPESTATUS[0]}${PIPESTATUS[1]}${PIPESTATUS[2]}" = 000 ]' _ {} \
    || { echo "per-contig restriction failed"; exit 1; }
  bcf concat --naive -Oz -o "$OUTDIR/$SAMPLE.targeted.vcf.gz" $(for c in "${CHROMS[@]}"; do echo "$PARTS/$c.vcf.gz"; done) \
    || { echo "concat failed"; exit 1; }
  read -r HOMREF_ROWS ALT_FILLED SKIP_MULTI DROP_HOMREF DROP_DUP < <(cat "$PARTS"/*.fill | awk '{a+=$1; b+=$2; c+=$3; d+=$4; e+=$5} END{print a+0, b+0, c+0, d+0, e+0}')
  rm -rf "$PARTS"
  printf '{"vcf_input_alt_fill": {"hom_ref_target_rows": %s, "alt_filled": %s, "skipped_multirecord": %s, "dropped_homref_dup": %s, "dropped_dup": %s, "multi_alt_policy": "first table ALT"}}\n' \
    "$HOMREF_ROWS" "$ALT_FILLED" "$SKIP_MULTI" "$DROP_HOMREF" "$DROP_DUP" > "$OUTDIR/vcf_input.json"
  echo "[vcf] alt_filled=$ALT_FILLED of $HOMREF_ROWS; skipped_multirecord=$SKIP_MULTI; dropped_homref_dup=$DROP_HOMREF; dropped_dup=$DROP_DUP"
  bcf index -t -f "$OUTDIR/$SAMPLE.targeted.vcf.gz"
  T_RESTRICT=$((SECONDS - t))
  EXPR=$(python3 bin/vcf_input.py restrict-expr)
  bcf view --threads "$THREADS" -t "$EXPR" -Oz -o "$SCORE_IN.vcf.gz" "$OUTDIR/$SAMPLE.targeted.vcf.gz"
  bcf index -t -f "$SCORE_IN.vcf.gz"
  OUT_ROWS=$(bcf index -n "$SCORE_IN.vcf.gz")
  AF=$(bcf query -f '%CHROM\t%POS[\t%GT]\n' "$SCORE_IN.vcf.gz" | python3 bin/vcf_input.py afrac "$TARGETS")
  FRAC=${AF##*fraction=}
  echo "[vcf] input rows=$IN_ROWS targets=$NTARGETS restricted rows=$OUT_ROWS $AF restrict=${T_RESTRICT}s total=$((SECONDS - T0))s"
  echo "[vcf] sex=${SEX:-unknown}; sex chromosomes (chrX/chrY/chrM) EXCLUDED from scoring: pgsc_calc imports VCFs without sex, so chrX PGS loci are dropped (README: 'chrX loci lost')"
  ALLSITES=$(awk -v f="$FRAC" 'BEGIN{print (f >= 0.5) ? "yes" : "no"}')
  [ "$ALLSITES" = no ] && echo "[vcf] WARNING: input is not an all-sites VCF (explicit genotype at $FRAC of target positions); match rate will be low" >&2
  if [ -n "$CHECK" ]; then
    echo "check-input: index=$IDX input_records=$IN_ROWS target_positions=$NTARGETS restricted_records=$OUT_ROWS alt_filled=$ALT_FILLED hom_ref_target_rows=$HOMREF_ROWS skipped_multirecord=$SKIP_MULTI dropped_homref_dup=$DROP_HOMREF dropped_dup=$DROP_DUP duplicate_allele_ids=$(bcf query -f '%CHROM:%POS:%REF:%ALT\n' "$OUTDIR/$SAMPLE.targeted.vcf.gz" | uniq -d | wc -l) explicit_genotype_fraction=$FRAC all_sites=$ALLSITES sex=${SEX:-unknown} restricted_size=$(du -h "$SCORE_IN.vcf.gz" | cut -f1) restrict_wall=${T_RESTRICT}s"
    [ "$ALLSITES" = no ] && { echo "CHECK_INPUT_FAILED: input is not an all-sites VCF (explicit genotype at $FRAC of $NTARGETS target positions, need >= 0.5)"; exit 1; }
    rm -rf "$OUTDIR/check_input"; mkdir -p "$OUTDIR/check_input"
    docker run --rm -v /data:/data -w "$PWD" "$PLINK2_IMG" plink2 --vcf "$SCORE_IN.vcf.gz" \
      --allow-extra-chr --chr 1-22,X,Y,XY --max-alleles 2 --set-all-var-ids '@:#:$r:$a' \
      --new-id-max-allele-len 100 missing --make-pgen --out "$OUTDIR/check_input/$SAMPLE" > "$OUTDIR/check_input/plink2.log" 2>&1 \
      && echo "CHECK_INPUT_OK: plink2 imported $SCORE_IN.vcf.gz ($(grep -vc '^#' "$OUTDIR/check_input/$SAMPLE.pvar" 2>/dev/null) variants, sex=${SEX:-unknown})" \
      || { echo "CHECK_INPUT_FAILED: see $OUTDIR/check_input/plink2.log"; tail -3 "$OUTDIR/check_input/plink2.log"; exit 1; }
    exit 0
  fi
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

if [ -n "$FASTPATH" ]; then
  # --- fastpath: same outputs as the pgsc_calc block below, cached sample-independent work ------
  [ -n "$PANEL" ] || { echo "--fastpath needs --panel (ancestry calibration; the cache is built for it)"; exit 2; }
  [ -n "$SCOREFILE_CACHE" ] || { echo "--fastpath needs --scorefile-cache covering every --pgs id"; exit 2; }
  [ -f "$SCORE_IN.vcf.gz" ] || { echo "--fastpath needs the prepped VCF $SCORE_IN.vcf.gz"; exit 2; }
  rm -rf "$OUTDIR/score"
  FP_ARGS=(--sample "$SAMPLE" --prepped "$SCORE_IN.vcf.gz" --outdir "$OUTDIR" --pgs "$PGS"
           --scorefile-cache "$SCOREFILE_CACHE" --panel "$PANEL" --ref "$REF")
  # resources: POLY_FASTPATH_THREADS (default NPROC, else 8) -> --threads/--score-threads; POLY_FASTPATH_MEM_GB (default 100)
  FP_ARGS+=(--threads "${POLY_FASTPATH_THREADS:-${NPROC:-8}}" --score-threads "${POLY_FASTPATH_THREADS:-${NPROC:-8}}"
            --mem-gb "${POLY_FASTPATH_MEM_GB:-100}")
  [ -n "$FASTPATH_CACHE" ] && FP_ARGS+=(--cache "$(readlink -f "$FASTPATH_CACHE")")
  : > "$OUTDIR/score.log"
  python3 bin/fastpath.py score "${FP_ARGS[@]}" >> "$OUTDIR/score.log" 2>&1 \
    || { echo "fastpath scoring failed (no fallback to pgsc_calc) -- see $OUTDIR/score.log"; tail -5 "$OUTDIR/score.log"; exit 1; }
else
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
    run_nf "$OUTDIR/score.log" run pgscatalog/pgsc_calc -profile docker \
      --input "$OUTDIR/samplesheet.csv" --target_build GRCh38 --scorefile "$bsf/*.txt.gz" --min_overlap 0.1 "${ANC[@]}" \
      -c conf/rootless.config --outdir "$bo" -work-dir "$WORKDIR" -resume \
      || { echo "scoring bin $k failed/timed out -- see $OUTDIR/score.log"; exit 1; }
    BIN_OUT+=("$bo")
  done
  mkdir -p "$OUTDIR/score"
  python3 bin/merge_scores.py "$OUTDIR/score" "$SAMPLE" "${BIN_OUT[@]}" >> "$OUTDIR/score.log" 2>&1
else
  : > "$OUTDIR/score.log"
  run_nf "$OUTDIR/score.log" run pgscatalog/pgsc_calc -profile docker \
    --input "$OUTDIR/samplesheet.csv" --target_build GRCh38 "${SCORE_SEL[@]}" --min_overlap 0.1 "${ANC[@]}" \
    -c conf/rootless.config --outdir "$OUTDIR/score" -work-dir "$WORKDIR" "${RESUME[@]}" \
    || { echo "scoring failed/timed out -- see $OUTDIR/score.log"; exit 1; }
fi

fi

[ -n "$SEX" ] && echo "$SEX" > "$OUTDIR/score/sample_sex.txt"
# raw-VCF ALT-fill counters -> score/vcf_input.json; grade_pgs puts them in provenance.json parameters
[ -n "$VCF" ] && [ -f "$OUTDIR/vcf_input.json" ] && cp -f "$OUTDIR/vcf_input.json" "$OUTDIR/score/vcf_input.json"
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
