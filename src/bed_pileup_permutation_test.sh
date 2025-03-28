#!/bin/bash

# bed_pileup_permutation_test.sh
# Previously we have used an arbitrary threshold of n >= 3 to consider a recurrently amplified locus "interesting".*
# Here we define a statistical threshold for enriched loci as the 99.9% confidence interval at above which regions
# would be unlikely to be recurrently amplified assuming a uniform distribution of ecDNA across the mappable genome.

# Usage: 
# ped_pileup_permutation_test.sh -i [directory_containing_bed_files] -g hg38.genome -b blacklist.bed

# Assumes conda is installed, you have a conda environment containing bedtools
conda activate bedtools

## Parse command line args
input_dir=""
genome=""
blacklist=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_dir|-i) input_dir="$2"; shift ;;
        --genome|-g) genome="$2"; shift ;;
        --blacklist|-b) blacklist="$2"; shift ;;
        --output_dir|-o) output_dir="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done
# Check if required arguments are provided
if [[ -z $input_dir || -z "$genome" || -z "$blacklist" || -z "$output_dir" ]]; then
    echo "Usage: $0 --input_dir <dir> --genome <genome> --blacklist <file> --output_dir <dir>"
    exit 1
fi

## Permutation test
mkdir -p $output_dir
for i in {1..100}; do
    # perform 1 permutation: shuffle input regions, and get genome coverage.
    cat "$input_dir"/*.bed | \
    bedtools shuffle -excl $blacklist -g $genome -i - | \
    bedtools sort -i - | \
    bedtools genomecov -bga -g $genome -i - | \
    bedtools subtract -a - -b $blacklist  > \
    $output_dir/${i}.bdg
done
