# Scripts for visualizing the number and distribution of ecDNA amplifications across the genome.

Source code for **Fig. 2, Suppl. Fig. 3 , Suppl. Tbl. 5** 

## Contents
`bed-pileup.ipynb` Plot distribution of ecDNA across the genome in .bedgraph format (**Fig. 2a**).  
`check-oncogenes.ipynb` Stringent filtering of oncogenes annotated in AmpliconClassifier results. (**Suppl. Tbl. 5**).  
`copy-number-comparisons.ipynb` Boxplots of ecDNA copy number by amplification type (**Fig. 2c-d**) and oncogene evidence (**Suppl. Fig. 3**).
`generate-bed-table.ipynb` Copy ecDNA genomic regions to tab-separated format (**Suppl. Tbl. 13**)
`gene-statistics.ipynb` Counting oncogenes most frequently amplified on ecDNA. 
`recurrent-amps.ipynb` More statistics counting ecDNA with and without oncogenes. 
`recurrent-amp-significance-threshold` Statistical basis for 'recurrent amplification' (**Suppl. Fig. 2, Suppl. Note 1**).
`recurrent-amps-plot.ipynb` Plot of recurrently amplified oncogenes by cancer type (**Fig. 2b**).
    `genes.txt` Config file for `reccurent-amps-plot`.
`symlink-bed-files.ipynb` Create symlinks of ecDNA regions as in `bed-pileup`, but organized by tumor and amplicon type.

## Requirements
See individual notebooks.

## Dependency tree
```
bed-pileup
    run_bed_pileup_permutation_test.sh
        recurrent-amp-significance-threshold
            recurrent-amps-plot
check-oncogenes
    data_imports.generate_gene_table()
        copy-number-comparisons
        gene-statistics
        recurrent-amps
symlink-bed-files
generate-bed-table
```