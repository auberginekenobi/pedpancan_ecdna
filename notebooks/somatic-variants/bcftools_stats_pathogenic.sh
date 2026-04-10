#!/bin/bash
set -euo pipefail

# Requirements: run bcftools_filter_pathogenic.sh first.
# Using bcftools, iterate through all CBTN .vcf.gz files,
# Filter for HotSpotAllele==1 (pathogenic)
# and write the bcftools stats outputs to a new file.

INPUT_DIR="./data/filter/pathogenic"
OUTPUT_DIR="./data/stats/pathogenic"

mkdir -p $OUTPUT_DIR

for file in $INPUT_DIR/*.vcf.gz; do
	OUT_BASE=$(basename $file .vcf.gz)
	OUT_FILE=$OUTPUT_DIR/${OUT_BASE}.stats.txt
	echo "Processing $OUT_FILE ..."
	bcftools stats -f PASS,. -i INFO/HotSpotAllele==1 $file > $OUT_FILE
done

