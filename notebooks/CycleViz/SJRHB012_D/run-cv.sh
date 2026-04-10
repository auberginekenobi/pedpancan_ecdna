#!/bin/bash
# Developed for CycleViz 0.2.1

source ~/.zshrc
conda activate cycleviz

SAMPLE=SJRHB012_D
CV_SRC=$HOME/software/CycleViz

# amplicon1
# amplicon1 has lots of discordant pairs that don't obviously form a cycle or correspond with CN changes, but a simple cycle accounts for most 
# of the amplification of ETV6.

CYCLES=../../../data/source/AmpliconArchitect/stjude/solid_paired/SJRHB012_D/SJRHB012_D_amplicon1_cycles.txt
GRAPH=../../../data/source/AmpliconArchitect/stjude/solid_paired/SJRHB012_D/SJRHB012_D_amplicon1_graph.txt
$CV_SRC/convert_cycles_file.py -c $CYCLES -g $GRAPH

NUMS=( 1 3 4 )
for i in "${NUMS[@]}"; do
$CV_SRC/CycleViz.py --cycles_file $CYCLES --cycle $i --graph $GRAPH \
        --ref GRCh38 --figure_size_style small \
	--gene_highlight_list MDM2 CDKN1B
done
