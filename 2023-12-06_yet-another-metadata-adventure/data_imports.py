#! python3
# Library of functions to load data from various places.
# 
# Owen Chapman

import pandas as pd
import pathlib
import os

## Functions to load metadata from the CAVATICA API. 
## See also 2023-11-27_cavatica-api/cavatica-api.ipynb
def clean_x01_biosample_metadata(df):
    '''
    Clean known errors in the x01 metadata.
    TODO: could clean the race and ethnicity fields but we don't use these.
    '''
    df = df.replace({
        'Tumor Descriptor':{
            "initial CNS Tumor": "Initial CNS Tumor",
            "Not Applicable":"Unavailable",
            "Diagnosis": "Initial CNS Tumor"
        },
        'gender':{
            "Not Reported":None
        }
    })
    # remove suffix from pnoc sample ids
    df.sample_id = df.sample_id.map(lambda x: '-'.join(x.split('-')[:2]) if x.startswith("7316-") else x)
    return df 
def import_x01_biosample_metadata(path="../2023-11-27_cavatica-api/out/X01-biosample-metadata.tsv"):
    path = pathlib.Path(path)
    df = pd.read_csv(path, sep='\t',index_col=0)
    df = clean_x01_biosample_metadata(df)
    df["cohort"]="PBTA-X01"
    return df
def import_x00_biosample_metadata(path="../2023-11-27_cavatica-api/out/X00-biosample-metadata.tsv"):
    df = import_x01_biosample_metadata(path)
    df["cohort"]="PBTA-X00"
    return df
def import_pnoc_biosample_metadata(path="../2023-11-27_cavatica-api/out/PNOC-biosample-metadata.tsv"):
    df = import_x01_biosample_metadata(path)
    df["cohort"]="PNOC"
    return df

## Function to load metadata from the AmpliconClassifier results
def get_pedpancan_biosamples_from_AC(include_x01=False,path='../data/local/AmpliconClassifier/pedpancan_summary_map.txt'):
    path = pathlib.Path(path)
    df = pd.read_csv(path, sep='\t', header=None, index_col=0, names = ["biosample","file"])
    if not include_x01:
        df['firstletter']=df.index.map(lambda x: x[0])
        df = df[df.firstletter != 'P']
    return df.index

## Function to compile CAVATICA metatdata for all CBTN samples in our cohort.
def import_cbtn_biosample_metadata(include_X01=False):
    if include_X01:
        df = pd.concat([import_x00_biosample_metadata(),import_x01_biosample_metadata(),import_pnoc_biosample_metadata()])
    else:
        df = pd.concat([import_x00_biosample_metadata(),import_pnoc_biosample_metadata()])
    cohort = get_pedpancan_biosamples_from_AC()
    df = df[df.index.isin(cohort)]
    return df

## Functions to open & preprocess opentarget histology data.
def clean_opentarget_histologies_files(df):
    cohort = import_cbtn_biosample_metadata()
    df = df[df.sample_id.isin(cohort.sample_id)]
    df = df[df.sample_type == 'Tumor'] # Drop normals
    df = df[df.composition != 'Derived Cell Line'] # Drop cell lines
    df = df[df.experimental_strategy != "Targeted Sequencing"] # these metadata are very different
    df = df.drop(["RNA_library","seq_center","pathology_free_text_diagnosis","gtex_group","gtex_subgroup","normal_fraction",
                  "cell_line_composition","cell_line_passage","tumor_fraction_RFpurify_ABSOLUTE",
                  "tumor_fraction_RFpurify_ESTIMATE","tumor_fraction_LUMP","dkfz_v12_methylation_mgmt_status",
                  "dkfz_v12_methylation_mgmt_estimated","integrated_diagnosis",
                  "tumor_fraction","tumor_ploidy","cohort"],axis=1) # drop columns we know we don't want
    df = df.replace({
        'composition':{
            "Not Available": None,
        },
        'extent_of_tumor_resection':{
            "Not Reported":None,
            'Unavailable':None,
            'Not Applicable':None
        },
    })
    
    # correct known errors
    df.loc["BS_K14VJ1E3","age_at_diagnosis_days"] = 2778
    df = df.drop(["BS_03G6PJKJ","BS_HJJPT3NR","BS_15R0SQRN","BS_GGDMSB26","BS_VQGR0D61"]) # lots of biosamples for the same tumor at the same timepoint
    
    # Propagate metadata from the same sample_id
    g = df.groupby('sample_id')
    df = []
    for name, group in g:
        columns = [col for col in group.columns if col not in ['sample_id','aliquot_id','experimental_strategy']]
        for column in columns:
            unique_values = group[column].dropna().unique()
            if len(unique_values) == 0:
                continue
            elif len(unique_values) == 1:
                non_na_value = unique_values[0]
                group[column].fillna(non_na_value, inplace=True)
            else:
                print(f"Warning: The column '{column}' for sample {name} differs between CAVATICA and opentarget annotations.")
        group=group.sort_values('experimental_strategy')
        df.append(group)
    df = pd.concat(df)
    
    # Add entries missing a KF biospecimen ID, but with a matching external biosample id.
    missing_bs = (cohort[~cohort.index.isin(df.index)]["sample_id"]).sort_values()
    print(f"{len(missing_bs)} KF biospecimens missing from the opentarget histologies table...")
    missing_bs = missing_bs[missing_bs.isin(df.sample_id)]
    print(f"found {len(missing_bs)} matching external sample IDs, adding to table...")
    newdf = []
    for biospecimen, sample in missing_bs.items():
        newentry=df[df.sample_id==sample].iloc[0]
        newentry.name = biospecimen
        newdf.append(newentry)
    newdf = pd.DataFrame(newdf)
    df = pd.concat([df,newdf])
    
    # Subset our cohort
    df = df[df.index.isin(cohort.index)]
    return df
def import_opentarget_histologies_files(path='/Users/ochapman/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/projects/2023-pedpancan/data/metadata/histologies.tsv'):
    path = pathlib.Path(path)
    df = pd.read_csv(path,sep='\t',index_col=0)
    df = clean_opentarget_histologies_files(df)
    return df

## Import Sunita's data
def import_sunita_classifications(path='../data/combinedamplicons.xlsx'):
    df = pd.read_excel(path)
    df = df[["sample_ID","cancer_type"]]
    df = df.drop_duplicates()
    df = df.set_index("sample_ID")
    return df

## Integrate all CBTN data available
def generate_cbtn_biosample_table():
    '''
    '''
    df = pd.DataFrame(index=get_pedpancan_biosamples_from_AC())
    cavatica_data = import_cbtn_biosample_metadata()
    cavatica_data = cavatica_data[['gender','Kids First Participant ID','cohort','disease_type','sample_id','Tumor Descriptor','primary_site','age_at_diagnosis']]
    df = pd.merge(left=df,how='inner',right=cavatica_data,left_index=True,right_index=True)
    opentarget_data = import_opentarget_histologies_files()
    #
    opentarget_data = opentarget_data.drop(["sample_id","composition","Kids_First_Participant_ID","experimental_strategy","sample_type",
                                           "germline_sex_estimate","race","ethnicity","molecular_subtype_methyl"],axis=1)
    #opentarget_data = opentarget_data[[]]
    # don't include race, ethnicity
    df = pd.merge(left=df,how='left',right=opentarget_data,left_index=True,right_index=True)
    
    return df

## SJ data
def import_sj_sample_info(path="/Users/ochapman/projects/pedpancan_ecdna/2022-02-23_sj_samples/SAMPLE_INFO_2022-03-02.tsv"):
    path = pathlib.Path(path)
    df = pd.read_csv(path,sep='\t',index_col="sample_name")
    return df