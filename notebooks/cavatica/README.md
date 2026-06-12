# Notebooks for accessing CBTN data using the sevenbridges API.

**These notebooks download intermedate results files using a private API key and are therefore not reproducible by the general user.**  

## Contents
`cavatica-api.ipynb` pull AmpliconSuite (aka AmpliconArchitect, AA) results.  
`download-mosdepth-results.ipynb` pull mosdepth results for CBTN samples.  
`download-variants.ipynb` pull consensus simple (SNV and indel) variants.  
`variant-manifests.ipynb` match AA results to extant variant calls.

## Configuration

Install the sevenbridges API:
```
conda env create -f sevenbridges.yml
NAME="sevenbridges"
conda activate $NAME
python -m ipykernel install --user --name '${NAME}' --display-name '${NAME}'
conda deactivate
```
Configure your credentials file at `~/.sevenbridges/credentials`:
```
[cavatica]
api_endpoint = https://cavatica-api.sbgenomics.com/v2
auth_token = [[ copy your authentication token from https://cavatica.sbgenomics.com/developer/token ]]
```
