#!/usr/bin/env bash
# poly-suite genotype-prep stage (chunk-parallel).
# Force-genotype a BAM at PGS scoring loci so hom-ref sites become dosage-0 (0/0)
# calls with the ALT allele retained (-A) — the fix for the variant-only-VCF
# coverage cap. Covered non-variant target -> 0/0; uncovered -> ./. (honest missing).
#
# usage: genotype_prep.sh <bam> <ref.fa> <out.vcf.gz> <harmonized_scorefile.txt.gz>...
# env:   NPROC (default 24) parallel chunks; BCFTOOLS_IMG to override the container.
set -euo pipefail

BAM=$(readlink -f "$1"); REF=$(readlink -f "$2"); OUT=$(readlink -f "$3"); shift 3
SCOREFILES=("$@")
IMG=${BCFTOOLS_IMG:-quay.io/biocontainers/bcftools:1.21--h8b25389_0}
NPROC=${NPROC:-24}
HERE=$(cd "$(dirname "$0")" && pwd)
OUTDIR=$(dirname "$OUT"); BASE=$(basename "${OUT%.vcf.gz}")
CHUNK="$OUTDIR/${BASE}_chunks"; rm -rf "$CHUNK"; mkdir -p "$CHUNK"
ALL="$OUTDIR/${BASE}.targets.tsv"

echo "[genotype-prep] building targets from ${#SCOREFILES[@]} scorefile(s)..."
python3 "$HERE/scoring_targets.py" "$ALL" "$REF" "${SCOREFILES[@]}"
NSITES=$(wc -l < "$ALL")
echo "[genotype-prep] $NSITES targets -> $NPROC parallel chunks"

# split into NPROC coord-ordered chunks (split -n l/N preserves line order, and
# ALL is already sorted, so chunk k..k+1 stays genome-sorted -> concat needs no re-sort)
split -n "l/$NPROC" -d -a 3 "$ALL" "$CHUNK/part_"

# force-genotype one chunk (alleles file for `call -C alleles`, regions for mpileup -R)
force_chunk() {
  local t=$1
  cut -f1,2 "$t" > "$t.regions"
  bgzip -f "$t";         tabix -s1 -b2 -e2 -f "$t.gz"
  bgzip -f "$t.regions"; tabix -s1 -b2 -e2 -f "$t.regions.gz"
  docker run --rm -v /data:/data -w "$PWD" "$IMG" bash -c "
    set -euo pipefail
    bcftools mpileup -f '$REF' -R '$t.regions.gz' -a FORMAT/DP,FORMAT/AD -Ou '$BAM' |
    bcftools call -m -A -C alleles -T '$t.gz' -Oz -o '$t.vcf.gz'
  "
}
export -f force_chunk; export REF BAM IMG

echo "[genotype-prep] force-genotyping..."
printf '%s\n' "$CHUNK"/part_* | xargs -P "$NPROC" -I{} bash -c 'force_chunk "$@"' _ {}

# concat chunks in order (each is coord-sorted; parts are in genome order)
echo "[genotype-prep] concatenating..."
mapfile -t PARTS < <(printf '%s\n' "$CHUNK"/part_*.vcf.gz | sort)
docker run --rm -v /data:/data -w "$PWD" "$IMG" \
  bcftools concat -Oz -o "$OUT" "${PARTS[@]}"
tabix -f -p vcf "$OUT"
rm -rf "$CHUNK"
echo "[genotype-prep] done -> $OUT"
echo "[genotype-prep] $(zcat "$OUT" | grep -vc '^#') genotype records emitted"
