#!/bin/bash
# Developed for CycleViz 0.2.1

source ~/.zshrc
conda activate cycleviz

SAMPLE=SJRHB012_S
CV_SRC=$HOME/software/CycleViz

# amplicon1
# amplicon1 has lots of discordant pairs that don't obviously form a cycle or correspond with CN changes, but a simple cycle accounts for most 
# of the amplification of ETV6.

CYCLES=../../../data/source/AmpliconArchitect/stjude/aa_workflow/SJRHB012_S/SJRHB012_S_amplicon1_cycles.txt
GRAPH=../../../data/source/AmpliconArchitect/stjude/aa_workflow/SJRHB012_S/SJRHB012_S_amplicon1_graph.txt

#$CV_SRC/convert_cycles_file.py -c $CYCLES -g $GRAPH

# Cycle 1 is represented in S.F. 9
NUMS=( 1 2 5 )
for i in "${NUMS[@]}"; do
$CV_SRC/CycleViz.py --cycles_file $CYCLES --cycle $i --graph $GRAPH \
        --ref GRCh38 --figure_size_style small
done

# amplicon2
# 2 cycles sharing most of their sequence account for all copy number amp of MDM2. The locations where the sequences differ each account for about half
# of the total amplification.

CYCLES=../../../data/source/AmpliconArchitect/stjude/aa_workflow/SJRHB012_S/SJRHB012_S_amplicon2_cycles.txt
GRAPH=../../../data/source/AmpliconArchitect/stjude/aa_workflow/SJRHB012_S/SJRHB012_S_amplicon2_graph.txt

#$CV_SRC/convert_cycles_file.py -c $CYCLES -g $GRAPH

NUMS=( 1 2 )
for i in "${NUMS[@]}"; do
$CV_SRC/CycleViz.py --cycles_file $CYCLES --cycle $i --graph $GRAPH \
	--ref GRCh38 --figure_size_style small \
	--gene_highlight_list MDM2
done

