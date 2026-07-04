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

SAMPLE= PGS= OUTDIR= VCF= BAM= PANEL= SEX= BOOT= DRY= REUSE_PREP= WORK_CACHE= SCOREFILE_CACHE=
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
ANC=(); [ -n "$PANEL" ] && ANC=(--run_ancestry "$PANEL")
FMT=vcf   # score-input format; reuse-prep upgrades to pfile when a .pgen exists (skips VCF->pgen)

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
  BOOT=${BOOT:-$VCF}
  [ -n "$BOOT" ] || { echo "BAM input needs --bootstrap-vcf (a VCF to trigger scorefile harmonization)"; exit 2; }
  plan "step 1: harmonize scorefiles (bootstrap $BOOT)"
  plan "step 2: genotype-prep BAM at scoring loci -> $OUTDIR/${SAMPLE}.prepped.vcf.gz"
  SCORE_IN="$OUTDIR/${SAMPLE}.prepped"
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
  printf 'sampleset,path_prefix,chrom,format\n%s,%s,,vcf\n' "$SAMPLE" "${BOOT%.vcf.gz}" > "$OUTDIR/hm.samplesheet.csv"
  rm -rf "$OUTDIR/hm" "$OUTDIR/work_hm"
  NXF_ANSI_LOG=false "$NF" run pgscatalog/pgsc_calc -profile docker \
    --input "$OUTDIR/hm.samplesheet.csv" --target_build GRCh38 --pgs_id "$PGS" --min_overlap 0.1 \
    -c conf/rootless.config --outdir "$OUTDIR/hm" -work-dir "$OUTDIR/work_hm" > "$OUTDIR/harmonize.log" 2>&1
  mapfile -t SF < <(python3 -c "import glob,os;s={};[s.setdefault(os.path.basename(p),p) for p in glob.glob('$OUTDIR/work_hm/**/normalised_*_hmPOS_GRCh38.txt.gz',recursive=True)];print(chr(10).join(sorted(s.values())))")
  [ "${#SF[@]}" -ge 1 ] || { echo "no harmonized scorefiles — check $OUTDIR/harmonize.log"; exit 1; }
  NPROC=${NPROC:-24} bin/genotype_prep.sh "$BAM" "$REF" "$OUTDIR/${SAMPLE}.prepped.vcf.gz" "${SF[@]}"
  [ -z "$SEX" ] && SEX=$(python3 bin/infer_sex.py "$BAM" | sed -n 's/.*sex=\([a-z]*\).*/\1/p')
fi

printf 'sampleset,path_prefix,chrom,format\n%s,%s,,%s\n' "$SAMPLE" "$SCORE_IN" "$FMT" > "$OUTDIR/samplesheet.csv"
rm -rf "$OUTDIR/score"

# --scorefile-cache: use locally-cached harmonized scoring files (--scorefile) instead of
# --pgs_id, skipping the ~40-58s DOWNLOAD_SCOREFILES step. Falls back to --pgs_id on any miss.
SCORE_SEL=(--pgs_id "$PGS")
if [ -n "$SCOREFILE_CACHE" ]; then
  SFDIR="$OUTDIR/scorefiles"; rm -rf "$SFDIR"; mkdir -p "$SFDIR"; miss=
  for id in ${PGS//,/ }; do
    f="$SCOREFILE_CACHE/${id}_hmPOS_GRCh38.txt.gz"
    if [ -f "$f" ]; then ln -sf "$(readlink -f "$f")" "$SFDIR/"; else miss="$miss $id"; fi
  done
  if [ -z "$miss" ]; then SCORE_SEL=(--scorefile "$SFDIR/*.txt.gz")
  else echo "scorefile-cache miss:$miss -> using --pgs_id" >&2; fi
fi
# --work-cache: persistent shared work dir + -resume, so the 16GB panel extraction
# (EXTRACT_DATABASE, keyed on the panel only) is reused across runs instead of re-unpacked.
WORKDIR="$OUTDIR/work"; RESUME=()
if [ -n "$WORK_CACHE" ]; then WORKDIR=$(readlink -f "$WORK_CACHE"); mkdir -p "$WORKDIR"; RESUME=(-resume)
else rm -rf "$OUTDIR/work"; fi
NXF_ANSI_LOG=false "$NF" run pgscatalog/pgsc_calc -profile docker \
  --input "$OUTDIR/samplesheet.csv" --target_build GRCh38 "${SCORE_SEL[@]}" --min_overlap 0.1 "${ANC[@]}" \
  -c conf/rootless.config --outdir "$OUTDIR/score" -work-dir "$WORKDIR" "${RESUME[@]}" > "$OUTDIR/score.log" 2>&1

[ -n "$SEX" ] && echo "$SEX" > "$OUTDIR/score/sample_sex.txt"
python3 bin/grade_pgs.py "$OUTDIR/score"
echo "RUN_DONE -> $OUTDIR/score  (pgs_scores.tsv/.json + provenance.json + report.html)"
