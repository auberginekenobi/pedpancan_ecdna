#!/bin/bash
# Developed for CycleViz 0.2.1

source ~/.zshrc
conda activate cycleviz

SAMPLE=BS_M4E4H6NG
CV_SRC=/Users/ochapman/software/CycleViz

CYCLES="/Users/ochapman/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/projects/2023-pedpancan/data/AmpliconArchitect/cbtn/hgg/PT_Y76A7PBA/BS_M4E4H6NG/BS_M4E4H6NG_amplicon1_cycles.txt"
GRAPH="/Users/ochapman/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/projects/2023-pedpancan/data/AmpliconArchitect/cbtn/hgg/PT_Y76A7PBA/BS_M4E4H6NG/BS_M4E4H6NG_amplicon1_graph.txt"
CYCLE=2

$CV_SRC/CycleViz.py --cycles_file $CYCLES --cycle $CYCLE --graph $GRAPH \
        --ref GRCh38 --figure_size_style small --outname $SAMPLE \
	--gene_highlight_list MYCL AKT3

