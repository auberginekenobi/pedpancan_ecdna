#! python3
# Download the SJ RNA-seq feature counts that the 2026-09-01 data request lists for S.T.2 biosamples, but that we don't yet hold.
# Usage, from notebooks/sjcloud in the dxpy env, logged in to DNAnexus (dx login) or with DX_AUTH_TOKEN set:
#   python download-sj-feature-counts.py             # list the files that would be downloaded
#   python download-sj-feature-counts.py --download  # download them
#
# Owen Chapman

import argparse
import os
import sys
import pandas as pd
import dxpy

sys.path.append('../../src')
import data_imports

PROJECT_ID = 'project-JBPg48j90qV46xFZ4K07k745' # ncu-pedpancan
SAMPLE_INFO = '../../data/source/sjcloud/SAMPLE_INFO_20260901.txt'
FC_DIR = '../../data/source/gene-expression/sjcloud/FEATURE_COUNTS'

def get_missing_feature_counts(sample_info, fc_dir):
    '''
    Return the feature count files in sample_info for S.T.2 biosamples which are not already in fc_dir.
    Files listed under more than one file_id (distinct RNA libraries with the same biosample name) are reported and skipped,
    since they would download to the same file name.
    '''
    biosamples = data_imports.import_biosamples()
    df = pd.read_csv(sample_info, sep='\t', low_memory=False)
    df = df[(df.sequencing_type == 'RNA-Seq') & (df.file_type == 'FEATURE_COUNTS') & df.sample_name.isin(biosamples.index)]
    df = df.drop_duplicates('file_id') # some files are listed once per dataset
    df['file_name'] = df.file_path.str.split('/').str[-1]
    df = df[~df.file_name.isin(os.listdir(fc_dir))]
    ambiguous = df.file_name.duplicated(keep=False)
    if ambiguous.any():
        print(f'Skipping {ambiguous.sum()} files with the same name but different file_id; choose one manually:')
        print(df[ambiguous][['sample_name','file_id','file_size','sj_datasets','attr_library_selection_protocol','attr_read_length']].to_string(index=False))
    return df[~ambiguous][['sample_name','file_id','file_name','file_size']]

def download_feature_counts(files, fc_dir, project_id):
    '''
    Download files to fc_dir. Write to a .part file and check its size before renaming, so an interrupted download is not
    mistaken for a complete one.
    '''
    for row in files.itertuples():
        part = os.path.join(fc_dir, row.file_name + '.part')
        dxpy.download_dxfile(row.file_id, part, project=project_id)
        if os.path.getsize(part) != row.file_size:
            raise IOError(f'{row.file_name}: downloaded {os.path.getsize(part)} bytes, expected {row.file_size}')
        os.rename(part, os.path.join(fc_dir, row.file_name))
        print(f'downloaded {row.file_name}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', action='store_true', help='download the files; otherwise only list them')
    args = parser.parse_args()

    files = get_missing_feature_counts(SAMPLE_INFO, FC_DIR)
    print(f'{len(files)} feature count files to download:')
    print(files.to_string(index=False))
    if args.download:
        if os.getenv('DX_AUTH_TOKEN'):
            dxpy.set_security_context({"auth_token_type": "Bearer", "auth_token": os.getenv('DX_AUTH_TOKEN')})
        download_feature_counts(files, FC_DIR, PROJECT_ID)