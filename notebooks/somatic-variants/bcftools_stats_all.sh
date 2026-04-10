#!/bin/bash
set -euo pipefail

# Using bcftools, iterate through all .vcf.gz files,
# and write the bcftools stats outputs to a new file.

DATA_DIR="../../data/variants"
CBTN_DIR="$DATA_DIR/cbtn_somatic"
SJ_DIR="$DATA_DIR/sj_somatic"
OUTPUT_DIR="./data/stats/all"

mkdir -p $OUTPUT_DIR

for file in $CBTN_DIR/*.vcf.gz $SJ_DIR/*.vcf.gz; do
	OUT_FILE=$(basename $file .vcf.gz)
	OUT_FILE=$OUTPUT_DIR/${OUT_FILE}.stats.txt
	echo "Processing $OUT_FILE ..."
	bcftools stats -f PASS,. $file > $OUT_FILE
done

