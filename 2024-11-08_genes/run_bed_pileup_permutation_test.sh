#!/bin/bash

# Define inputs
AA_DATA_REPO=$HOME/anno/AmpliconArchitect/GRCh38
genome=$HOME/anno/genomes/human.hg38.noalt.genome
blacklist=$AA_DATA_REPO/GRCh38_merged_centromeres_conserved_sorted.bed
bed_dir=$(pwd)/bed_symlinks/ecDNA_all

# Run
perm_test="../src/bed_pileup_permutation_test.sh"
$perm_test -i $bed_dir -g $genome -b $blacklist -o "tmp/ecDNA"

bed_dir=$(pwd)/bed_symlinks/intrachromosomal_all

$perm_test -i $bed_dir -g $genome -b $blacklist -o "tmp/intrachromosomal"