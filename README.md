# Source code for analyses of ecDNA in pediatric cancer (pedpancan_ecDNA)
Contains working code, scripts and data for the Chavez Lab pedpancan project.  
Read our [preprint on medRxiv](http://doi.org/10.1101/2025.07.22.24308163)!  

Tested on:
- Apple M2 Pro chip and 16Gb RAM running macOS Sonoma 14.5.
- Apple M3 Pro chip and 18Gb RAM running macOS Sequoia 15.6.1.

## Installation

*Install time:* 20 minutes.

### `conda` installation

(Almost) all code is in jupyter notebook format. Packages and dependencies are installed using `conda`, and are specified in `.yml` files in the `env` directory.  For instructions on how to set up `jupyter` on your workstation, see steps 1-5 of [Setting up your workstation](https://github.com/auberginekenobi/protocols/tree/main/0_Setting_up_your_workstation). Dependencies are indicated in the first cell of each notebook. To install a conda environment from a .yml file, run
```
## Create a new environment and install all packages
NAME="myenvironment"
conda env create -f ${NAME}.yml
## OR ##
## If you're on a Mac with Apple silicon, R packages need to be installed using intel architecture:
CONDA_SUBDIR=osx-64 conda env create -f ${NAME}.yml
conda activate ${NAME}
conda config --env --set subdir osx-64


## Link the environment to your base jupyter installation
# R environments:
conda activate ${NAME}
Rscript -e "IRkernel::installspec(name = '${NAME}', displayname = '${NAME}')"
conda deactivate

# python environments:
conda activate $NAME
python -m ipykernel install --user --name ${NAME} --display-name ${NAME}
conda deactivate
```

### Other dependencies
- [CycleViz](https://github.com/AmpliconSuite/CycleViz) v0.2.1  
- [oscutils](https://github.com/auberginekenobi/oscutils) @adc5b13 or later 

Source directories for these dependencies should be downloaded to `~/software`.

### Required source data

To reproduce figures and statistical tests, users should download the following:
- Suppl. Tbls. for this publication
- [AmpliconSuite results for this publication](https://ampliconrepository.org/project/69c59d0cd69472646656266f)
- [Suppl. Tbls. from Chapman *et al*, 2023](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01551-3/MediaObjects/41588_2023_1551_MOESM4_ESM.xlsx)
- [Suppl. Tbls. 2, 5 from Corbett *et al*, 2025](https://www.nature.com/articles/s41467-025-65190-4#Sec30)
- [AmpliconArchitect Data Repo](https://github.com/AmpliconSuite/AmpliconArchitect#setting-up-the-aa-data-repo)

The `data` directory of this repository should be organized as follows:[^1]
```
data
├── Supplementary Tables.xlsx
├── annot
│   ├── AmpliconArchitect
│   │   └── GRCh38
│   │   │   └── [.. AA_DATA_REPO files..]
│   └── human.hg38.noalt.genome
├── external
│   ├── Chapman2023
│   │   └── 41588_2023_1551_MOESM4_ESM.xlsx
│   └── Corbett2025
│       ├── 41467_2025_65190_MOESM4_ESM.xlsx
│       └── 41467_2025_65190_MOESM7_ESM.xlsx
└── source
    ├── AmpliconArchitect
    │   └── [... AmpliconArchitect output files ...]
    └── AmpliconClassifier
        └── [... AmpliconClassifier output files ...]
```
## Usage

This repository contains source code necessary to reproduce figures and statistical tests from the **Supplementary Tables** and other publicly available data. 

### Reproducing figures and statistical tests

*Running time: 1 hour.*

Having successfully installed dependencies and downloaded required data in the section above, the following notebooks may be run using `jupyter lab`:
```
notebooks
├── CycleViz # Suppl. Figs. 4c, 9
│   ├── BS_M4E4H6NG
│   │   └── run-cv.sh 
│   ├── SJHGG052_A
│   │   └── run-cv.sh 
│   ├── SJRHB012_D
│   │   └── run-cv.sh
│   ├── SJRHB012_S
│   │   └── run-cv.sh 
│   ├── SJRHB031519_D1
│   │   └── run-cv.sh 
│   └── cycle_lengths.ipynb 
│
├── genes-genomic-regions # Run in this order
│   ├── bed-pileup.ipynb # Fig. 2a
│   ├── run_bed_pileup_permutation_test.sh # Suppl. Fig. 2
│   ├── recurrent-amp-significance-threshold.ipynb # Suppl. Fig. 1, Suppl. Note 1
│   ├── recurrent-amps-plot.ipynb # Fig. 2b
│   ├── copy-number-comparisons.ipynb # Fig. 2c-d, Suppl. Fig. 4
│   ├── gene-statistics.ipynb # statistics
│   ├── recurrent-amps.ipynb # statistics
│   └── generate-bed-table.ipynb # Suppl. Tbl. 13
│
├── germline-variants
│   └── germline-associations.ipynb # Suppl. Note 5
│
├── longitudinal-analyses
│   └── longitudinal-samples.ipynb # Fig. 5
│
├── metadata-analyses
│   ├── figure1btables.Rmd # Fig. 1b
│   ├── medullo-diff.ipynb # Suppl. Note 1
│   ├── notes-by-tumor-type.ipynb # various Results
│   └── summary_statistics.ipynb # various Results
│
├── misc-qc
│   └── age_regressions.ipynb # Suppl. Fig. 1
│
└── survival
    └── survival.ipynb # Fig. 4, Suppl. Figs. 6-7
```

Code to generate the **Suppl. Tbls.** is also included in this repository, but requires source data which are access-controlled, under embargo, licensed, or otherwise not freely available. Notebooks not on the above list require input data that are not currently publicly available. If running these analyses is of interest, please contact the corresponding author.

## Dataset
Most source data for this publication are organized in the Supplementary Tables. Code for generating and reading the Suppl. Tbls. is in [data_imports.py](src/data_imports.py). 

To read the Suppl. Tbls.:
```
# import data_imports.py
import sys
sys.path.append('../src')
from data_imports import *

patients = import_patients()
biosamples = import_biosamples()
amplicons = import_amplicons()
genes = import_genes()
```

To generate the Suppl. Tbls. from source data:
```
biosamples = generate_biosamples_table()
patients = generate_patient_table(biosamples)
amplicons = generate_amplicon_table(biosamples)
genes = generate_gene_table(biosamples)
```
To generate the Suppl. Tables, the following source files are required:
- data/source/AmpliconClassifier/pedpancan_summary_map.txt # list of all biosamples analyzed. Generated by AmpliconClassifier/ampclasslib/make_input.py.
- data/source/AmpliconClassifier/pedpancan_amplicon_classification_profiles.tsv # Amplicon classifications. Generated by AmpliconClassifier/amplicon_classifier.py.
- data/source/AmpliconClassifier/pedpancan_gene_list.tsv # Amplified genes. Generated by ibid.
- data/Supplementary Tables.xlsx sheet '9. Tumor ontology' # Compiled by the authors from [St Jude](https://permalinks.stjude.cloud/permalinks/st-jude-cloud-disease-ontology) and [DKFZ](https://www.molecularneuropathology.org/mnp/classifiers/11) ontologies.
- data/source/sjcloud/SAMPLE_INFO_SJ00.txt # File metadata generated by the St. Jude Cloud upon file provision.
- data/source/opentarget/histologies.tsv # File metadata from the OpenPBTA project ([source](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/analyses/molecular-subtyping-integrate/results/histologies.tsv)).
- data/source/cavatica/X01-biosample-metadata.tsv # File metadata for CBTN dataset. Compiled using the CAVATICA API. See [cavatica-api.ipynb](notebooks/cavatica/cavatica-api.ipynb).
- data/source/cavatica/X00-biosample-metadata.tsv # Ibid.
- data/source/cavatica/PNOC-biosample-metadata.tsv # Ibid.
- data/oncogenes/oncogene_blacklist.txt # List of genes with insufficient oncogenic evidence, generated by [check-oncogenes.ipynb](notebooks/genes-genomic-regions/check-oncogenes.ipynb).
- data/external/Dubois2022/NIHMS1907773-supplement-Supplemental_tables_1-6.xlsx # Suppl. Tbls. from Dubois *et al*, 2022, used to annotate some shared samples.

## Contributions
`main` branch contains working code for our most recent manuscript iteration (4/2026), and `dev` branch contains latest incremental code updates. To contribute, please follow this workflow:
- Checkout a new branch from `dev`. Name it something descriptive.
- Develop on your branch.
- **Do not push data, large files, or images, or notebooks with the same embedded, to git.** If you need to share these files, use the project OneDrive. The reason for this is that git keeps permanent copies of all files pushed to the repo. Thus, even if you later delete a file, it still takes up space in the git history. This would be merely inconvenient in the case of large files, images, etc. but would become a major issue if sensitive data (access tokens, patient data) were added. Below are some conventions; ask Owen if you have any questions.
  - I usually put input data files in `./data`. `./data` is in the `.gitignore` to prevent these files being added to the repo.
  - Output files (figures, intermediate data files, etc.) often go in `./myanalysisfolder/out`. This folder is likewise in the `.gitignore`.
  - Jupyter notebooks display cell outputs and other information which are useful during development but not source code. Thus, before committing a Jupyter notebook I usually do 'Kernel' > 'Restart kernel and clear outputs of all cells' before saving and committing.
- When you are ready to merge, make a pull request to `dev`.
- Do a code review with Owen.
- Complete the pull request, delete the old branch.

[^1]: Currently, AmpliconSuite results downloaded from AmpliconRepository need to be moved and renamed to match expected paths in notebooks. It's next on the todo list.
