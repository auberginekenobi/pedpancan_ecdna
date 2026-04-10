#!/bin/bash
# Developed for CycleViz 0.2.1

source ~/.zshrc
conda activate cycleviz

SAMPLE=E-1336-21
CV_SRC=/Users/ochapman/software/CycleViz

CYCLES="../../manuscript/collaborators/Konstantin/E-1336-21/E-1336-21_AA_results/E-1336-21_amplicon1_cycles.txt"
GRAPH="../../manuscript/collaborators/Konstantin/E-1336-21/E-1336-21_AA_results/E-1336-21_amplicon1_graph.txt"
CYCLE=1

$CV_SRC/CycleViz.py --cycles_file $CYCLES --cycle $CYCLE --graph $GRAPH \
        --ref GRCh38 --figure_size_style small --outname $SAMPLE \
	--gene_highlight_list MYCN

