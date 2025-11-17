#!/bin/bash

# Using bcftools, iterate through all .vcf.gz files,
# and write the bcftools stats outputs to a new file.

INPUT_DIR="./data/filter/sj_likely_pathogenic"
OUTPUT_DIR="./data/sj_stats/likely_pathogenic"

for file in $INPUT_DIR/*.vcf.gz; do
	OUT_FILE=$(basename $file .consensus_somatic.norm.annot.public.filtered.vcf.gz)
	OUT_FILE=$OUTPUT_DIR/${OUT_FILE}.stats.txt
	echo "Processing $OUT_FILE ..."
	bcftools stats -f PASS,. $file > $OUT_FILE
done
