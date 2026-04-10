#!/bin/bash
## Run permutation test to evaluate statistical likelihood of independent overlapping ecDNA regions.
## See also src/bed_pileup_permutation_test.sh

# Define inputs
PROJ_ROOT=../..
AA_DATA_REPO=$PROJ_ROOT/anno/AmpliconArchitect/GRCh38
genome=$PROJ_ROOT/anno/genomes/human.hg38.noalt.genome
blacklist=$AA_DATA_REPO/GRCh38_merged_centromeres_conserved_sorted.bed
bed_dir=$(pwd)/bed_symlinks/ecDNA_all

# Run
perm_test="../../src/bed_pileup_permutation_test.sh"
$perm_test -i $bed_dir -g $genome -b $blacklist -o "tmp/ecDNA"

bed_dir=$(pwd)/bed_symlinks/intrachromosomal_all

$perm_test -i $bed_dir -g $genome -b $blacklist -o "tmp/intrachromosomal"