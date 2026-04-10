#!/bin/bash


source /home/ochapman/.bashrc
mamba activate CycleViz-mamba


SAMPLE=SJRHB012_D
CV_SRC=$lukas_home/bin/CycleViz

# amplicon1
# amplicon1 has lots of discordant pairs that don't obviously form a cycle or correspond with CN changes, but a simple cycle accounts for most 
# of the amplification of ETV6.

CYCLES=/expanse/lustre/projects/csd677/collab/projects/pedpancan/AmpliconArchitect/stjude/solid_paired/SJRHB012_D/SJRHB012_D_amplicon1_cycles.txt
GRAPH=/expanse/lustre/projects/csd677/collab/projects/pedpancan/AmpliconArchitect/stjude/solid_paired/SJRHB012_D/SJRHB012_D_amplicon1_graph.txt
$CV_SRC/convert_cycles_file.py -c $CYCLES -g $GRAPH

NUMS=( 1 3 4 )
for i in "${NUMS[@]}"; do
$CV_SRC/CycleViz.py --cycles_file ${SAMPLE}_amplicon1_BPG_converted_cycles.txt --cycle $i \
        --ref GRCh38 --figure_size_style small \
	--gene_highlight_list MDM2 CDKN1B
done
