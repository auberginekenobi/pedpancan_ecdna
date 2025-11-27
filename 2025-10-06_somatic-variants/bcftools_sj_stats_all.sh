#!/bin/bash

# Using bcftools, iterate through all .vcf.gz files,
# and write the bcftools stats outputs to a new file.

INPUT_DIR="./data/sj_somatic"
OUTPUT_DIR="./data/sj_stats/all"

for file in $INPUT_DIR/*.vcf.gz; do
	OUT_FILE=$(basename $file .consensus_sj_somatic.norm.annot.public.vcf.gz)
	OUT_FILE=$OUTPUT_DIR/${OUT_FILE}.stats.txt
	echo "Processing $OUT_FILE ..."
	bcftools stats -f PASS,. $file > $OUT_FILE
done