#!/usr/bin/env bash
# Score the auto-selected best-evidence PGS set (bin/select_pgs.py ->
# results/pgs_catalog_meta.json) end-to-end on the HG001 BAM:
#   pass 1  harmonize the selected scorefiles (throwaway score on benchmark VCF)
#   prep    genotype-prep the BAM at the union of the selected loci
#   pass 2  re-score on the prepped full-genotype VCF  (good coverage)
# Grade separately:  python3 bin/grade_pgs.py results/selected
set -euo pipefail
cd /data/alvin/poly-suite

BAM=/data/alvin/ref/GIAB/HG001.bwa.sortdup.bqsr.bam
REF=/data/alvin/ref/GRCh38/hg38.canonical.fa
IDS=$(python3 -c "import json;print(','.join(sorted(json.load(open('results/pgs_catalog_meta.json')))))")
echo "[selected] auto-selected set: $IDS"

echo "[selected] pass 1 — harmonize scorefiles"
rm -rf results/selected_hm work_selected_hm
NXF_ANSI_LOG=false /home/alvin/bin/nextflow run pgscatalog/pgsc_calc -profile docker \
  --input data/samplesheet.csv --target_build GRCh38 \
  --pgs_id "$IDS" --min_overlap 0.1 \
  -c conf/rootless.config --outdir results/selected_hm -work-dir work_selected_hm \
  > selected_pass1.log 2>&1

mapfile -t SF < <(python3 - <<'PY'
import glob, os
seen = {}
for p in glob.glob("work_selected_hm/**/normalised_*_hmPOS_GRCh38.txt.gz", recursive=True):
    seen.setdefault(os.path.basename(p), p)   # one per PGS id
for b in sorted(seen):
    print(seen[b])
PY
)
echo "[selected] harmonized scorefiles: ${#SF[@]}"
[ "${#SF[@]}" -ge 1 ] || { echo "NO_SCOREFILES — pass 1 failed"; exit 1; }

echo "[selected] genotype-prep at union of selected loci (~15M)"
NPROC=32 bin/genotype_prep.sh "$BAM" "$REF" results/prep/HG001_selected.vcf.gz "${SF[@]}"

echo "[selected] pass 2 — re-score on prepped VCF"
printf 'sampleset,path_prefix,chrom,format\nHG001sel,/data/alvin/poly-suite/results/prep/HG001_selected,,vcf\n' \
  > data/samplesheet_selected.csv
rm -rf results/selected work_selected
NXF_ANSI_LOG=false /home/alvin/bin/nextflow run pgscatalog/pgsc_calc -profile docker \
  --input data/samplesheet_selected.csv --target_build GRCh38 \
  --pgs_id "$IDS" --min_overlap 0.1 \
  -c conf/rootless.config --outdir results/selected -work-dir work_selected \
  > selected_pass2.log 2>&1

python3 bin/infer_sex.py "$BAM" results/selected >/dev/null 2>&1 || true
echo "SELECTED_DONE exit=0"
