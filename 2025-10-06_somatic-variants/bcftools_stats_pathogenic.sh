#!/bin/bash

INPUT_DIR="./data/somatic"
OUTPUT_DIR="./data/stats/pathogenic"

for file in $INPUT_DIR/*.vcf.gz; do
	OUT_FILE=$(basename $file .consensus_somatic.norm.annot.public.vcf.gz)
	OUT_FILE=$OUTPUT_DIR/${OUT_FILE}.stats.vcf
	echo "Processing $OUT_FILE ..."
	bcftools stats -f PASS,. -i INFO/HotSpotAllele==1 $file > $OUT_FILE
done

