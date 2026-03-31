#!/bin/bash
set -euo pipefail

INPUT_DIR="../data/variants/cbtn_somatic"
OUTPUT_DIR="./data/filter/pathogenic"

mkdir -p $OUTPUT_DIR

for file in $INPUT_DIR/*.vcf.gz; do
	OUT_BASE=$(basename $file .consensus_somatic.norm.annot.public.vcf.gz)
	OUT_FILE=$OUTPUT_DIR/${OUT_BASE}.p.vcf.gz
	echo "Processing $OUT_FILE ..."
	bcftools filter -i 'INFO/HotSpotAllele==1 && (FILTER="PASS" || FILTER=".")' -Oz -o $OUT_FILE $file
    bcftools index -t $OUT_FILE
	#break
done


