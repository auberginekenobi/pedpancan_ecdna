# Notebooks for accessing data using various APIs.

## Sevenbridges API

Setup:
- Install the sevenbridges API:
```
conda env create -f sevenbridges.yml
NAME="sevenbridges"
conda activate $NAME
python -m ipykernel install --user --name '${NAME}' --display-name '${NAME}'
conda deactivate
```
- Configure your credentials file at `~/.sevenbridges/credentials`:
```
[cavatica]
api_endpoint = https://cavatica-api.sbgenomics.com/v2
auth_token = [[ copy your authentication token from https://cavatica.sbgenomics.com/developer/token ]]
```
