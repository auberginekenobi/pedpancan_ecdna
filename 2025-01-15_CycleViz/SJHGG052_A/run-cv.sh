#!/bin/bash
# Developed for CycleViz 0.2.1

source ~/.zshrc
conda activate cycleviz

SAMPLE=SJHGG052_A
CV_SRC=/Users/ochapman/software/CycleViz

CYCLES="/Users/ochapman/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/projects/2023-pedpancan/data/AmpliconArchitect/stjude/results_batch_20240912/SJHGG052_A/SJHGG052_A_amplicon1_cycles.txt"
GRAPH="/Users/ochapman/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/projects/2023-pedpancan/data/AmpliconArchitect/stjude/results_batch_20240912/SJHGG052_A/SJHGG052_A_amplicon1_graph.txt"
CYCLE=7

$CV_SRC/CycleViz.py --cycles_file $CYCLES --cycle $CYCLE --graph $GRAPH \
        --ref GRCh38 --figure_size_style small --outname $SAMPLE \
	--gene_highlight_list AKT3

