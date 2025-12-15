#!/bin/bash
set -euo pipefail

# Wrapper script runs bcftools_filter_likely_pathogenic.sh

DATA_DIR="../data/variants"
CBTN_DIR="$DATA_DIR/cbtn_somatic"
SJ_DIR="$DATA_DIR/sj_somatic"
OUTPUT_DIR="./data/filter/likely_pathogenic"
GENES="./cancer_genes.txt"
FILTER_SCRIPT="./bcftools_filter_likely_pathogenic.sh"

mkdir -p $OUTPUT_DIR

for file in $CBTN_DIR/*.vcf.gz $SJ_DIR/*.vcf.gz; do
    base=$(basename "$file" .vcf.gz)
    OUT_FILE="$OUTPUT_DIR/${base}.lp.vcf.gz"
    echo "Processing $file  ..."
    "$FILTER_SCRIPT" -i "$file" -g "$GENES" -o $OUT_FILE
done

echo "All samples filtered. Output in $OUTPUT_DIR"