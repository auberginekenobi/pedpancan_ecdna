#!/bin/bash
set -euo pipefail
INPUT_DIR="./data/somatic"
OUTPUT_DIR="./data/filter/likely_pathogenic"
mkdir -p $OUTPUT_DIR
GENES="/Users/gotoriki/pedpancan_ecdna/2025-10-06_somatic-variants/cancer_genes.txt"
FILTER_SCRIPT="/Users/gotoriki/pedpancan_ecdna/2025-10-06_somatic-variants/bcftools_filter_path_permissive.sh"

for file in "$INPUT_DIR"/*.vcf.gz; do
    base=$(basename "$file" .vcf.gz)
    OUT_FILE="$OUTPUT_DIR/${base}.filtered.vcf.gz"
    echo "Processing $file -> $OUT_FILE ..."
    "$FILTER_SCRIPT" -i "$file" -g "$GENES" -o "$OUT_FILE"
done

echo "All samples filtered. Output in $OUTPUT_DIR"