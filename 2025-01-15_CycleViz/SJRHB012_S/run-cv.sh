#!/bin/bash


source /home/ochapman/.bashrc
mamba activate CycleViz-mamba


SAMPLE=SJRHB012_S
CV_SRC=$lukas_home/bin/CycleViz

# amplicon1
# amplicon1 has lots of discordant pairs that don't obviously form a cycle or correspond with CN changes, but a simple cycle accounts for most 
# of the amplification of ETV6.

CYCLES=/expanse/lustre/projects/csd677/collab/projects/pedpancan/AmpliconArchitect/stjude/aa_workflow/SJRHB012_S/SJRHB012_S_amplicon1_cycles.txt
GRAPH=/expanse/lustre/projects/csd677/collab/projects/pedpancan/AmpliconArchitect/stjude/aa_workflow/SJRHB012_S/SJRHB012_S_amplicon1_graph.txt

#$CV_SRC/convert_cycles_file.py -c $CYCLES -g $GRAPH
# Modified the output to use cycle 14; not sure why this isn't included...

$CV_SRC/CycleViz.py --cycles_file ${SAMPLE}_amplicon1_BPG_converted_cycles.txt --cycle 14 \
        --ref GRCh38 --figure_size_style small


# amplicon2
# 2 cycles sharing most of their sequence account for all copy number amp of MDM2. The locations where the sequences differ each account for about half
# of the total amplification.

CYCLES=/expanse/lustre/projects/csd677/collab/projects/pedpancan/AmpliconArchitect/stjude/aa_workflow/SJRHB012_S/SJRHB012_S_amplicon2_cycles.txt
GRAPH=/expanse/lustre/projects/csd677/collab/projects/pedpancan/AmpliconArchitect/stjude/aa_workflow/SJRHB012_S/SJRHB012_S_amplicon2_graph.txt

$CV_SRC/convert_cycles_file.py -c $CYCLES -g $GRAPH

NUMS=( 1 2 )
for i in "${NUMS[@]}"; do
$CV_SRC/CycleViz.py --cycles_file ${SAMPLE}_amplicon2_BPG_converted_cycles.txt --cycle $i \
	--ref GRCh38 --figure_size_style small \
	--gene_highlight_list MDM2
done

