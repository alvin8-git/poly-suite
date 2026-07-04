#!/usr/bin/env bash
# Pre-convert a cached prepped VCF to plink2 pfile (.pgen/.pvar/.psam) ONCE, so
# reuse-prep runs point pgsc_calc at the pfile (format=pfile) and skip the slow
# VCF->pgen conversion (pgsc_calc's MAKE_COMPATIBLE) on every subsequent run.
# Variant IDs set to chr:pos:ref:alt to match pgsc_calc's matching scheme.
#
# usage: make_pgen.sh <prepped.vcf.gz>  ->  <prefix>.pgen/.pvar/.psam
set -euo pipefail
VCF=$(readlink -f "$1")
PREFIX="${VCF%.vcf.gz}"
IMG=${PLINK2_IMG:-ghcr.io/pgscatalog/plink2:2.00a5.10}

docker run --rm -v /data:/data -w "$PWD" "$IMG" \
  plink2 --vcf "$VCF" --make-pgen \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 \
    --out "$PREFIX" >/dev/null 2>&1
echo "pgen -> $PREFIX.{pgen,pvar,psam}  ($(grep -vc '^#' "$PREFIX.pvar" 2>/dev/null || echo '?') variants)"
