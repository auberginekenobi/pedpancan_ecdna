# pedpancan_ecdna
Contains working code, scripts and data for the Chavez Lab pedpancan project. Not yet public.  
Tested on an Apple M2 Pro chip and 16Gb RAM running macOS Sonoma 14.5.  

## Installation

(Almost) all code is in jupyter notebook format. Packages and dependencies are installed using `conda`.  For instructions on how to set up conda on your workstation, see [Setting up your workstation](https://github.com/chavez-lab/protocols/tree/main/Setting_up_your_workstation). Dependencies should be clearly indicated in the first cell of each notebook. To install a conda environment from a .yml file, run
```
## Create a new environment and install all packages
conda env create -f environment.yml

## If you're on a Mac with Apple silicon, edgeR and others need to be installed using intel architecture:
CONDA_SUBDIR=osx-64 conda env create -f differential-expression.yml
conda activate differential-expression
conda config --env --set subdir osx-64

## Link the environment to your base jupyter installation
# R environments:
conda activate myenvironment
Rscript -e "IRkernel::installspec(name = 'myenvironment', displayname = 'myenvironment')"
conda deactivate

# python environments:
conda activate myenvironment
python -m ipykernel install --user --name myenvironment --display-name "myenvironment"
conda deactivate
```

## TODO from S. Danovi
- Include multiome data analysis
- Is it ecDNA or just amplification?
