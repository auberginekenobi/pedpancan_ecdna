#!/bin/bash
set -euo pipefail

# Using bcftools, iterate through previously filtered .vcf.gz files,
# and write the bcftools stats outputs to a new file.

INPUT_DIR="./data/filter/likely_pathogenic"
OUTPUT_DIR="./data/stats/likely_pathogenic"

mkdir -p $OUTPUT_DIR

for file in $INPUT_DIR/*.vcf.gz; do
	OUT_FILE=$(basename $file .vcf.gz)
	OUT_FILE=$OUTPUT_DIR/${OUT_FILE}.stats.txt
	echo "Processing $OUT_FILE ..."
	bcftools stats -f PASS,. $file > $OUT_FILE
done
