# somatic SNV library, to be used with the various notebooks in the same directory.

import pandas as pd
import os

def parse_bcftools_stats(file):
    # read a bcftools stats output, and
    # return the total number of variants
    try:
        with open(file) as f:
            for line in f:
                if line.startswith('SN	0	number of records:'):
                    count=line.rstrip().split('\t')[-1]
                    count=int(count)
        return count
    except:
        print("could not read file:",file)
        raise

def get_variant_counts(path):
    # for all .stats.txt files in path,
    # create a dictionary of the format:
    # return {sample : count}
    d = {}
    for filename in os.listdir(path):
        full_path = os.path.join(path, filename)
        if full_path.endswith('.stats.txt'):
            count = parse_bcftools_stats(full_path)
            #print(count)
            filename = filename.split('.')[0]
            d[filename] = count
    return d

def get_variant_count_df():
    # return a dataframe of the format:
    # file_id, all_counts, likely_pathogenic_counts, pathogenic_counts

    columns = ['all','likely_pathogenic','pathogenic']
    for i in range(len(columns)):
        col = columns[i]
        counts = get_variant_counts(f'./data/stats/{col}')
        s = pd.Series(counts)
        s.name = f"{col}_counts"
        if i == 0:
            df = pd.DataFrame(s)
        else:
            df = df.join(s)
    return df

def read_cbtn_manifest(file='../../data/variants/metadata/manifest_20250910_143948.csv'):
    df = pd.read_csv(file)
    df=df[df.name.str.endswith('.gz')]
    df['index'] = df.name.map(lambda x: os.path.basename(x).split('.')[0])
    df = df.set_index('index')
    df = df.rename(columns={'Kids First Biospecimen ID':'biosample_id'})
    df = df['biosample_id'].copy()
    return df
def read_sj_manifest(file='../../data/variants/metadata/sj_somatic_vcfs.tsv'):
    df = pd.read_table(file)
    df = df[df['vcf_file_name'] != 'SJST030131_D3_G1.Somatic.vcf.gz'] # this file is corrupted
    df['index'] = df["vcf_file_name"].str.rsplit(".", n=-1).str[0]
    df = df.set_index('index')
    df = df.rename(columns={'sample_name':'biosample_id'})
    return df['biosample_id'].copy()
def read_manifest():
    d1 = read_cbtn_manifest()
    d2 = read_sj_manifest()
    df = pd.concat([d1,d2])
    check_file_bs_ids(df)
    return df

def reindex_counts(manifest,counts):
    df = counts.join(manifest,how='inner').set_index('biosample_id')
    return df

def check_file_bs_ids(example):
    # check no duplicates
    assert len(example.index) == len(example.index.unique())
    assert len(example.values) == len(example.unique())
    # check no missing values
    assert sum(example.index.isna()) == 0
    assert sum(example.isna()) == 0
    return

def deduplicate_snv_set(df):
    # deduplicate multiple biosamples with somatic snv data
    # df: dataframe with required columns ['patient_id','amplicon_class','amplicon_class', 'ecDNA_sequences_detected',
    # 'age_at_diagnosis','external_sample_id','file_name','ecDNA_sequences_detected','in_unique_tumor_set',
    # 'in_unique_patient_set']
    df = df.sort_values(by=['patient_id','amplicon_class','ecDNA_sequences_detected','age_at_diagnosis'],
                   ascending=[True,False,True,True])
    df["in_snv_set"]=~df.duplicated(subset=["patient_id"],keep='last')
    return df

def reannotate_snv_set(df):
    # add other useful annotations, remove the not useful ones
    df=df.copy()
    df["amplicon_class"] = df["amplicon_class"].replace({
        "intrachromosomal": "chromosomal",
    })
    df['amplified'] = df.amplicon_class.isin(['ecDNA','chromosomal'])
    df['ecDNA'] = df.amplicon_class == 'ecDNA'
    df = df.drop(columns=['external_sample_id','file_name','ecDNA_sequences_detected','in_unique_tumor_set','in_unique_patient_set'])
    #df = df[df.in_snv_set].copy()
    return df

def merge_ecDNA_counts(ecDNA_df,counts_df):
    df = ecDNA_df.join(counts_df,how='inner')
    df = deduplicate_snv_set(df)
    df = reannotate_snv_set(df)
    return df

def get_target_tumor_types(df,n=5):
    # find all cancer_types in df with at least of each amplicon_class and at least 5 samples. 
    stats = df.groupby("cancer_type").agg(
        count=("amplicon_class", "size"),
        classes=("amplicon_class", "nunique")
    )
    targets = stats.loc[
        (stats['count'] >= n) &
        (stats['classes'] == df['amplicon_class'].nunique())
    ]
    return targets.index.values