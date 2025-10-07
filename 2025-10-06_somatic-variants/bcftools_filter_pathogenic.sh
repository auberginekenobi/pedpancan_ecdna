#!/bin/bash

INPUT_DIR="./data/somatic"
OUTPUT_DIR="./data/filter/pathogenic"

for file in $INPUT_DIR/*.vcf.gz; do
	OUT_FILE=$(basename $file .consensus_somatic.norm.annot.public.vcf.gz)
	OUT_FILE=$OUTPUT_DIR/${OUT_FILE}.filter.vcf
	echo "Processing $OUT_FILE ..."
	bcftools filter -i 'INFO/HotSpotAllele==1 && (FILTER="PASS" || FILTER=".")' $file > $OUT_FILE
	#break
done


