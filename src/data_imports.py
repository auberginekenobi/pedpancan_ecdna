#! python3
# Library of functions to load metadata from various places, and merge into unified patient, biosample tables.
# Usage: generate_biosample_table()
# Requires: 
# 
# Owen Chapman

import pandas as pd
import numpy as np
import pathlib
import os
import warnings

## Function to load metadata from the AmpliconClassifier results
## Get this file from /expanse/lustre/projects/csd677/collab/projects/pedpancan/AmpliconClassifier/batch/inputs
def get_pedpancan_biosamples_from_AC(path='../data/source/AmpliconClassifier/pedpancan_summary_map.txt'):
    path = pathlib.Path(path)
    df = pd.read_csv(path, sep='\t', header=None, index_col=0, names = ["biosample","file"])
    df = df[~df.index.duplicated(keep='first')] # drop duplicate AA runs
    return df.index

## Functions to load metadata from the CAVATICA API. 
## See also 2023-11-27_cavatica-api/cavatica-api.ipynb
def import_x01_biosample_metadata(path="../data/source/cavatica/X01-biosample-metadata.tsv"):
    path = pathlib.Path(path)
    df = pd.read_csv(path, sep='\t',index_col=0)
    df["cohort"]="PBTA-X01"
    return df
def import_x00_biosample_metadata(path="../data/source/cavatica/X00-biosample-metadata.tsv"):
    df = import_x01_biosample_metadata(path)
    df["cohort"]="PBTA-X00"
    return df
def import_pnoc_biosample_metadata(path="../data/source/cavatica/PNOC-biosample-metadata.tsv"):
    df = import_x01_biosample_metadata(path)
    df["cohort"]="PNOC"
    return df
def clean_cavatica_biosample_metadata(df):
    '''
    Clean known errors in the x01 metadata, and unify ontologies.
    '''
    # drop duplicates
    df = df[~df.index.duplicated(keep='first')]
    # remove suffix from pnoc sample ids
    df.loc[:,"sample_id"] = df.sample_id.map(lambda x: '-'.join(x.split('-')[:2]) if x.startswith("7316-") else x)
    
    df = df.replace({
        'Tumor Descriptor':{
            "initial CNS Tumor": "Diagnosis",
            "Not Applicable":np.nan,
            "Unavailable":np.nan,
            "Initial CNS Tumor": "Diagnosis",
            "Progressive Disease Post-Mortem":"Progressive",
            "Deceased":"Autopsy",
            "Deceased - No Sample Collection":"Autopsy",
        },
        'gender':{
            "Not Reported":np.nan
        }
    })
    # Correct suspected errors
    df.loc["BS_6Z213H2V","Tumor Descriptor"] = "Progressive" # This tumor was resected 175 days after the initial tumor resection.
    df.loc["BS_1135HC0V","Tumor Descriptor"] = "Second Malignancy" # This dysplasia was diagnosed 574 days after the first tumor was resected, in a new location.
    df.loc["BS_ZS1QRMXS","Tumor Descriptor"] = "Progressive" # Tumor resected 128 days after previous resection.
    df.loc["BS_FVYBGMG1","Tumor Descriptor"] = "Progressive" # Tumor resected 107 days after previous resection.
    df.loc["BS_5J5VH3X0","Tumor Descriptor"] = "Progressive" # Biopsied 240 days after previous biopsy.
    df.loc["BS_E9M7TDB6","Tumor Descriptor"] = "Progressive" # Second resection in different location 112 days after previous resection.
    df.loc["BS_5XZP7F4Q","Tumor Descriptor"] = "Progressive" # Second resection 1845 days after previous
    df.loc["BS_EXTEGB51","Tumor Descriptor"] = "Progressive" # Third resection 1922 days after previous
    df.loc["BS_93BV8AY9","Tumor Descriptor"] = "Second Malignancy" # Second diagnosis 2975 days after initial.
    df.loc["BS_CRKBDAYZ","Tumor Descriptor"] = "Progressive" # Series of progressive diagnoses long after initial.
    df.loc["BS_B4DY7ET3","Tumor Descriptor"] = "Progressive" # Second resection 119 days after previous.
    return df 
## Function to compile CAVATICA metatdata for all CBTN samples in our cohort.
def import_cbtn_biosample_metadata(include_X01=True):
    if include_X01:
        df = pd.concat([import_x01_biosample_metadata(),import_x00_biosample_metadata(),import_pnoc_biosample_metadata()])
    else:
        df = pd.concat([import_x00_biosample_metadata(),import_pnoc_biosample_metadata()])
    cohort = get_pedpancan_biosamples_from_AC()
    df = df[df.index.isin(cohort)]
    df = clean_cavatica_biosample_metadata(df)
    return df

## Functions to open & preprocess opentarget histology data.
## Get histologies.tsv from https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/analyses/molecular-subtyping-integrate/results/histologies.tsv
def clean_opentarget_histologies_files(df,verbose=False):
    '''
    If verbose, include various tumor purity estimates.
    '''
    cohort = import_cbtn_biosample_metadata()
    df = df[df.sample_id.isin(cohort.sample_id)]
    df = df[df.sample_type == 'Tumor'] # Drop normals
    df = df[df.experimental_strategy != "Targeted Sequencing"] # these metadata are very different
    df = df.drop(["RNA_library","seq_center","pathology_free_text_diagnosis","gtex_group","gtex_subgroup","normal_fraction",
                  "cell_line_composition","cell_line_passage","dkfz_v12_methylation_mgmt_status",
                  "dkfz_v12_methylation_mgmt_estimated","integrated_diagnosis",
                  "tumor_ploidy","cohort"],axis=1) # drop columns we know we don't want
    if not verbose:
        df = df.drop([
            "tumor_fraction", # Theta2 purity estimates from WGS
            "tumor_fraction_RFpurify_ABSOLUTE", # RFpurify a random forest on 450k methylation, trying to replicate ABSOLUTE (SNP data)
            "tumor_fraction_RFpurify_ESTIMATE", # now trying to replicate ESTIMATE (RNA-seq or expression array)
            "tumor_fraction_LUMP", # 450k leukocyte-specific methylation sites, doi: 10.1038/ncomms9971
        ],axis=1)
    df = df.replace({
        'composition':{
            "Not Available": np.nan,
        },
        'extent_of_tumor_resection':{
            "Not Reported":np.nan,
            'Unavailable':np.nan,
            'Not Applicable':np.nan
        },
        'harmonized_diagnosis':{
            "Not Reported":np.nan
        },
        'tumor_descriptor':{
            "Not Applicable":np.nan,
            "Unavailable":np.nan,
            "Initial CNS Tumor": "Diagnosis",
            "Progressive Disease Post-Mortem":"Progressive",
            "Deceased":"Autopsy",
            "Deceased - No Sample Collection":"Autopsy",
        },
    })
    
    # correct known errors
    df.loc["BS_K14VJ1E3","age_at_diagnosis_days"] = 2778
    
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
    
    # Propagate metadata from the same sample_id
    g = df.groupby('sample_id')
    df = []
    for name, group in g:
        columns = [col for col in group.columns if col not in ['sample_id','aliquot_id','experimental_strategy','composition']]
        for column in columns:
            unique_values = group[column].dropna().unique()
            if len(unique_values) == 0:
                continue
            # if only 1 unique value for this column, propagate to all rows.
            elif len(unique_values) == 1:
                non_na_value = unique_values[0]
                group[column] = group[column].fillna(non_na_value)
            # For methylation columns, take the most confident methylation classifier score.
            elif column in ['dkfz_v11_methylation_subclass','dkfz_v12_methylation_subclass']:
                continue
            elif column == 'dkfz_v11_methylation_subclass_score':
                max_idx = group[column].idxmax()
                group[column] = group.loc[max_idx,column]
                group['dkfz_v11_methylation_subclass'] = group.loc[max_idx,'dkfz_v11_methylation_subclass']
            elif column == 'dkfz_v12_methylation_subclass_score':
                max_idx = group[column].idxmax()
                group[column] = group.loc[max_idx,column]
                group['dkfz_v12_methylation_subclass'] = group.loc[max_idx,'dkfz_v12_methylation_subclass']
            # if more than 1 unique value, throw a warning and do not change the table.
            else:
                warnings.warn(f"The column '{column}' for sample {name} differs between CAVATICA and opentarget annotations: {unique_values}")
        group=group.sort_values('experimental_strategy')
        df.append(group)
    df = pd.concat(df)
    
    # Subset our cohort
    df = df[df.index.isin(cohort.index)]
    return df
def import_opentarget_histologies_files(path='../data/cloud/opentarget/histologies.tsv',verbose=False):
    '''
    Get this file from /Users/ochapman/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/projects/2023-pedpancan/data/opentarget/histologies.tsv
    '''
    path = pathlib.Path(path)
    df = pd.read_csv(path,sep='\t',index_col=0,low_memory=False)
    df = clean_opentarget_histologies_files(df,verbose=verbose)
    return df

def propagate(df,dest,source,rename=False):
    '''
    Replace NA values in dest with those in source, then drop source and rename dest
    '''
    df[dest] = df[dest].fillna(df[source])
    df.drop(source, inplace=True, axis=1)
    if rename:
        df = df.rename(columns={dest:rename})
    return df
def consensus(df,dest,source,rename=False):
    '''
    Check that values in dest and source agree, if not then set to NA. Drop source and rename dest.
    '''
    df.loc[df[dest] != df[source], dest] = pd.NA
    df.drop(source, inplace=True, axis=1)
    if rename:
        df = df.rename(columns={dest:rename})
    return df    

## Integrate all CBTN data available
def generate_cbtn_biosample_table(verbose=0):
    '''
    Generate a metadata table of cbtn biosamples.
    verbose:
        0: most useful metadata only. These are included in the biosamples metadata table (S.T.2)
        1: includes some extra columns. Useful for generating the patient metadata table (S.T.1)
        2: includes lots of extra columns
    '''
    df = pd.DataFrame(index=get_pedpancan_biosamples_from_AC())
    cavatica_data = import_cbtn_biosample_metadata()
    df = pd.merge(left=df,how='inner',right=cavatica_data,left_index=True,right_index=True)
    verbosity = verbose != 0
    opentarget_data = import_opentarget_histologies_files(verbose=verbosity)
    df = pd.merge(left=df,how='left',right=opentarget_data,left_index=True,right_index=True,suffixes=(None,"_y"))

    # For CAVATICA annotations which are missing, propagate those from opentarget.
    df = propagate(df,"primary_site","primary_site_y")
    df = propagate(df,"age_at_diagnosis","age_at_diagnosis_days")
    df = propagate(df,"Tumor Descriptor","tumor_descriptor")
    df = consensus(df,"gender","reported_gender","sex")
    
    # For selected biosamples missing annotations, propagate from other biosample from same tumor.
    df.loc['BS_AH3RVK53'] = df.loc['BS_AH3RVK53'].fillna(df.loc['BS_G65EA38C'])
    df.loc['BS_KQPCYZ2K'] = df.loc['BS_KQPCYZ2K'].fillna(df.loc['BS_4DYW3T2A'])
    df.loc['BS_JEZBA2EW'] = df.loc['BS_JEZBA2EW'].fillna(df.loc['BS_JDZX545X'])
    df.loc['BS_XNYQS1WG'] = df.loc['BS_XNYQS1WG'].fillna(df.loc['BS_JDZX545X'])
    
    # Rename columns
    df.index.name = "biosample_id"
    df = df.rename(columns={
        'Kids First Participant ID':'patient_id',
        'Tumor Descriptor':'tumor_history',
        'case_id':'external_patient_id',
        'sample_id':'external_sample_id',
    })
    
    # Drop cell lines
    df = df[df.composition != 'Derived Cell Line']

    # drop columns
    if verbose < 2:
        df = df.drop(["race","ethnicity","external_patient_id","WGS_UUID","Kids First Biospecimen ID Normal",
                      "sample_id_y","composition","Kids_First_Participant_ID","experimental_strategy","sample_type",
                      "germline_sex_estimate","race_y","ethnicity_y","molecular_subtype_methyl","cohort_participant_id","Notes"
                     ],axis=1)  
    if verbose < 1:
        df = df.drop(["primary_site","pathology_diagnosis","OS_days","OS_status","EFS_days","age_last_update_days","aliquot_id",
                      "cancer_predispositions","CNS_region","age_at_chemo_start","age_at_radiation_start","cancer_group",
                      "age_at_event_days","clinical_status_at_event"
                     ],axis=1)
    
    return df

## SJ data

#def import_sj_sample_info(path="../data/local/sjcloud/SAMPLE_INFO_2022-03-02.tsv"):
def import_sj_sample_info(path="../data/cloud/sjcloud/SAMPLE_INFO_SJ00.txt"):
    path = pathlib.Path(path)
    df = pd.read_csv(path,sep='\t',index_col="sample_name")
    return df
def clean_sj_biosample_metadata(df):
    '''
    Clean known errors in the sj metadata, and unify ontologies, units etc.
    '''
    # drop duplicates
    columns = ['subject_name','sample_type','attr_age_at_diagnosis','attr_sex','sj_long_disease_name','sj_diseases','attr_oncotree_disease_code','sj_dataset_accessions']
    order = [True]*7+[False]
    df = df[(df.sequencing_type == 'WGS') & (df.file_type == 'BAM') & df.file_path.str.endswith('.bam')]
    df = df.sort_values(by=columns,ascending=order)
    df = df[~df.index.duplicated(keep='first')]
    
    # clean metadata classes
    df.sample_type = df.sample_type.map(str.title)

    df = df.replace({
        'attr_age_at_diagnosis':{
            "Not Available": np.nan
        },
        'attr_sex':{
            "Not Available":np.nan
        },
        'sample_type':{
            "Relapse":"Recurrence",
        }
    })
    # Convert age from years to days
    df['attr_age_at_diagnosis'] = (pd.to_numeric(df['attr_age_at_diagnosis'],errors='coerce')*365.25).round()

    # Rename columns
    df = df.rename(columns={
        'subject_name':'patient_id',
        'sample_type':'tumor_history',
        'attr_sex':'sex',
        'sj_dataset_accessions':'cohort',
        'attr_age_at_diagnosis':'age_at_diagnosis',
    })
    return df

def import_dubois_supplementary_data(path='../data/cloud/sjcloud/NIHMS1907773-supplement-Supplemental_tables_1-6.xlsx'):
    df = pd.read_excel(path,header=1)
    # drop unused columns
    df = df.drop(['Tumor_Sample_Barcode_Long','Autopsy','N_SNV','total_codingSNV','N_SV','Publication Alias','other_published_sample_ID','Histone']
                ,axis=1)
    # rename columns
    df = df.rename(columns={
        'Tumor_Sample_Barcode':'biosample_id',
        'Surviva_Status':'OS_status',
        'Overall_Survival_Months':'OS_months'
    })
    # subset SJ samples
    df = df[df.biosample_id.str.startswith('SJHGG')]
    df = df.set_index('biosample_id')

    # standardize terms
    df = df.replace({
        'OS_status':{
            "DECEASED": "Deceased",
            "LIVING":"Alive",
        }
    })
    return df

def generate_sj_biosample_table(verbose=0):
    '''
    Notes:
    sj_diseases != attr_oncotree_disease_code = sj_associated_diagnoses_disease_code
    attr_diagnosis != sj_long_disease_name != sj_associated_diagnoses
    '''
    df = pd.DataFrame(index=get_pedpancan_biosamples_from_AC())
    add = import_sj_sample_info()
    add = clean_sj_biosample_metadata(add)
    df = pd.merge(left=df,how='inner',right=add, left_index=True, right_index=True)
    
    df.index.name = "biosample_id"

    # drop columns
    if verbose < 2:
        df = df.drop(["file_path","file_id","sequencing_type","file_type","description","sj_embargo_date","attr_ethnicity","attr_race",
                      "sj_genome_build","sj_pipeline_name","sj_pipeline_version","attr_library_selection_protocol","attr_read_length",
                      "attr_sequencing_platform","attr_read_type","attr_tissue_preservative","attr_inferred_strandedness",
                      "attr_lab_strandedness","attr_germline_sample"
        ],axis=1)
    if verbose < 1:
        df = df.drop(["sj_pmid_accessions","sj_publication_titles","sj_pub_accessions","sj_datasets","sj_ega_accessions","attr_diagnosis",
                      "attr_diagnosis_group","attr_oncotree_disease_code","attr_subtype_biomarkers","sj_associated_diagnoses",
                      "sj_associated_diagnoses_disease_code"
        ],axis=1)

    # add annotations from Dubois et al 2021 if available
    try:
        dubois = import_dubois_supplementary_data()
        df['dubois'] = dubois.Histone_group
    except:
        # if we can't add annotations from Dubois et al 2021, don't add them.
        pass
    
    return df

# Function to create the cancer_type column based on priority
def get_subtype(row):
    priority_columns = ['molecular_subtype','dkfz_v12_methylation_subclass',
                    'dkfz_v11_methylation_subclass', 'harmonized_diagnosis', 'disease_type', "dubois", "sj_diseases"]  # Add other columns as needed
    for col in priority_columns:
        if col in ['dkfz_v12_methylation_subclass', 'dkfz_v11_methylation_subclass'] and pd.notnull(row[col]) and row[f"{col}_score"] > 0.9:
            if row[col].startswith("CONTR") or row[col].startswith("CTRL"):
                continue
            else:
                return col+','+row[col]
        elif col not in ['dkfz_v12_methylation_subclass', 'dkfz_v11_methylation_subclass'] and pd.notnull(row[col]):
            return col+','+row[col]
    return None
def unify_tumor_diagnoses(df, include_HM=False, path="../data/source/pedpancan_mapping.xlsx"):
    # Apply the function to create the cancer_subtype column
    path = pathlib.Path(path)
    mapping = pd.read_excel(path, 'mapping')
    mapping_dict = dict(zip(mapping['source_ontology']+','+mapping['source_class'], mapping['target_class']))
    submap_dict = dict(zip(mapping['source_ontology']+','+mapping['source_class'], mapping['target_subclass']))
    df['cancer_type'] = df.apply(get_subtype, axis=1)
    df['cancer_subclass'] = df['cancer_type'].map(submap_dict)
    df['cancer_type'] = df['cancer_type'].map(mapping_dict)
    # drop tumor type annotations now that we have a unified diagnosis.
    df = df.drop(["disease_type","dkfz_v11_methylation_subclass","dkfz_v11_methylation_subclass_score",
                  "dkfz_v12_methylation_subclass","dkfz_v12_methylation_subclass_score","molecular_subtype","harmonized_diagnosis",
                  "broad_histology","short_histology", "sj_long_disease_name", "sj_diseases", "dubois"
                 ],axis=1,errors='ignore')
    # Drop nontumor samples
    if include_HM:
        df = df[~df.cancer_type.isin(["NONTUMOR","UNLABELLED"])]
    else:
        df = df[~df.cancer_type.isin(["NONTUMOR","UNLABELLED","HM"])]
    return df

def clean_tumor_diagnoses(df):
    '''
    Correct known errors in tumor type annotations.
    '''
    df.loc['BS_AQMKA8NC','cancer_type']='ETMR' # Use diagnosis of primary.
    return df

## Annotate with ecDNA status
def annotate_with_ecDNA(df,path="../data/source/AmpliconClassifier/pedpancan_amplicon_classification_profiles.tsv"):
    '''
    Annotate biosamples with ecDNA status.
    Inputs:
        df: pd.DataFrame. Must be indexed by biosample.
        path: path to AmpliconClassifier results.
    '''
    # load AC results
    if path.endswith("Supplementary Tables.xlsx"):
        ac = pd.read_excel(path,sheet_name="4. Amplicons")
    else:
        ac = pd.read_csv(path,sep='\t')
    

    # Aggregate by biosample
    ac_agg = ac.groupby("sample_name").sum().ecDNA_amplicons
    df = df.join(ac_agg)
    df = df.rename(columns={"ecDNA_amplicons":"ecDNA_sequences_detected"})
    df["ecDNA_sequences_detected"] = df["ecDNA_sequences_detected"].fillna(0)
    return df

def amplicon_class_priority(df):
    '''
    Return the highest-priority amplicon class for a df. 
    Use with pd.groupby.apply. See annotate_amplicon_class.
    Inputs:
        df: pd.Dataframe, with sample_name and amplicon_decomposition_class columns. sample_name should have only 1 value.
    Returns:
        string: amplicon_class of highest priority for that sample_name.
    '''
    ec = df['ecDNA+'].values
    bfb = df['BFB+'].values
    classes = df.amplicon_decomposition_class.values
    if 'Positive' in ec:
        return 'ecDNA'
    elif 'Positive' in bfb:
        return 'intrachromosomal'
    elif 'Complex-non-cyclic' in classes:
        return 'intrachromosomal'
    elif 'Linear' in classes:
        return 'intrachromosomal'
    elif 'Complex non-cyclic' in classes: # these latter 2 are deprecated classes but present in some published datasets
        return 'intrachromosomal'
    elif 'Linear amplification' in classes:
        return 'intrachromosomal'
    else:
        return 'no amplification'

## Annotate with amplicon class (ecDNA, BFB, complex noncircular, linear, or none in descending priority order)
def annotate_amplicon_class(df,path="../data/source/AmpliconClassifier/pedpancan_amplicon_classification_profiles.tsv"):
    '''
    Annotate biosamples with amplicon class.
    Inputs:
        df: pd.Dataframe, indexed by biosample
        path: path to AmpliconClassifier results.
    '''
    # load AC results
    if path.endswith("Supplementary Tables.xlsx"):
        ac = pd.read_excel(path,sheet_name="4. Amplicons")
    else:
        ac = pd.read_csv(path,sep='\t')
        

    ac_agg = ac.groupby("sample_name")[['ecDNA+', 'BFB+', 'amplicon_decomposition_class']].apply(amplicon_class_priority)
    ac_agg.name = 'amplicon_class'
    df = df.join(ac_agg)
    df["amplicon_class"] = df["amplicon_class"].fillna('no amplification')
    return df

def annotate_duplicate_biosamples(df):
    '''
    Annotate duplicate biosamples by patient (same patient id) and tumor (same patient and tumor type).
    Priority given to ecDNA+ samples, then most recent.
    NB. Use unique_patient_set for survival, unique_tumor_set for figure 3.
    HACK HACK the unique_patient_set and unique_tumor_set identities are determined by the sort order. In particular, changing the sort
    order of the amplicon classes (or the string values used to encode them in amplicon_class_priority() or annotate_amplicon_class()
    will change the biosamples used in survival analysis, etc. 
    '''
    df = df.sort_values(by=['patient_id','amplicon_class','ecDNA_sequences_detected','age_at_diagnosis','external_sample_id'],
                       ascending=[True,False,True,True,False])
    df["in_unique_tumor_set"]=~df.duplicated(subset=["cancer_type","patient_id"],keep='last')
    df["in_unique_patient_set"]=~df.duplicated(subset=["patient_id"],keep='last')
    df = df.sort_values(by=['patient_id','age_at_diagnosis'],ascending=True)
    return df
    
def generate_biosample_table(include_HM=False,):
    df = pd.concat([generate_cbtn_biosample_table(),generate_sj_biosample_table()])
    df = unify_tumor_diagnoses(df,include_HM=include_HM)
    df = clean_tumor_diagnoses(df)
    df = annotate_with_ecDNA(df)
    df = annotate_amplicon_class(df)
    df = annotate_duplicate_biosamples(df)
    return df

## Generate Suppl. Table 1
###

def import_sj_survival_data(path="../data/cloud/sjcloud/SJ_SurvivalMaster.xlsx"):
    path = pathlib.Path(path)
    df = pd.read_excel(path,index_col=0)
    return df
def clean_sj_survival_data(df):
    df = df.dropna(subset=['Date of Primary Dx']).copy()
    df['tmp']=df['Date of Death'].fillna(df['Date of data collection'])
    df['OS_months'] = (df.tmp - df['Date of Primary Dx']).apply(lambda x:x.days * 12 / 365.25)
    df = df.rename(columns={
        'Survival Status':'OS_status'
    })
    df = df[['OS_status','OS_months']]
    df = df.replace({
        'OS_status':{
            "Expired": "Deceased",
        }
    })
    return df
def import_clean_cbtn_survival_data():
    df = generate_cbtn_biosample_table(verbose=1)
    df['OS_months']=df['OS_days']*12/365.25
    df = df[['OS_status','OS_months']]
    df = df.replace({
        'OS_status':{
            "DECEASED": "Deceased",
            "LIVING":"Alive",
        }
    })
    return df

def cat_sj_dubois_survival(sj,dubois = None):
    '''
    SJ df indexed by biosample id, with columns [OS_status 	OS_months]
    Add dubois rows if not already present.
    '''
    if dubois is None:
        dubois = import_dubois_supplementary_data()
    # get columns of interest
    dubois = dubois[["OS_months","OS_status"]]
    # subset only rows not already in sj
    dubois = dubois[~dubois.index.isin(sj.index)]
    return pd.concat([sj,dubois])
    
def generate_patient_table(biosamples_tbl=None):
    # Start with biosamples
    warnings.filterwarnings('ignore', '.*differs between CAVATICA and opentarget annotations.*')
    if biosamples_tbl is None:
        biosamples_tbl = generate_biosample_table()
    df = biosamples_tbl.copy()
    df = df[df.in_unique_patient_set == True]
    df = df[['sex','patient_id','age_at_diagnosis','cohort','cancer_type','cancer_subclass','amplicon_class']]
    # Add sj survival data
    surv = import_sj_survival_data()
    surv = clean_sj_survival_data(surv)
    surv = cat_sj_dubois_survival(surv)
    # Add cbtn survival data
    surv = pd.concat([surv,import_clean_cbtn_survival_data()])
    df = df.join(surv)
    df.set_index('patient_id',inplace=True)
    warnings.resetwarnings()
    return df

def generate_amplicon_table(biosamples_tbl=None,
                            path='../data/source/AmpliconClassifier/pedpancan_amplicon_classification_profiles.tsv'):
    # get biosample table if it wasn't provided
    if biosamples_tbl is None:
        biosamples_tbl = generate_biosample_table()
    # generate amplicon table
    df = pd.read_csv(path,sep='\t')
    df = df[df.sample_name.isin(biosamples_tbl.index)]
    # check for duplicates
    dups = df[df.duplicated(subset=['sample_name','amplicon_number'],keep=False)]
    if len(dups) > 0:
        warnings.warn(f'Duplicate amplicon table entries detected: \n {dups.to_string()}')
    return df

def generate_gene_table(biosamples_tbl=None,
                        path='../data/source/AmpliconClassifier/pedpancan_gene_list.tsv',
                        oncogene_blacklist_file='../data/oncogenes/oncogene_blacklist.txt'):
    # get biosample table if it wasn't provided
    if biosamples_tbl is None:
        biosamples_tbl = generate_biosample_table()
    # generate amplicon table
    df = pd.read_csv(path,sep='\t')
    df = df[df.sample_name.isin(biosamples_tbl.index)]
    # check for duplicates
    dups = df[df.duplicated(subset=['sample_name','amplicon_number','gene','truncated','feature'],keep=False)]
    if len(dups) > 0:
        warnings.warn(f'Duplicate gene table entries detected: \n {dups.to_string()}')
    # Edit oncogene annotations using the blacklist (see check-oncogenes.ipynb)
    with open(oncogene_blacklist_file, "r") as f:
        blacklist = set(line.strip() for line in f)
    df.loc[df['gene'].isin(blacklist),'is_canonical_oncogene'] = False
    return df

###########################################

## Supplementary Table imports
SUPPLEMENTARY_TABLES_PATH="../data/Supplementary Tables.xlsx"

def import_patients():
    return pd.read_excel(SUPPLEMENTARY_TABLES_PATH,sheet_name="1. Patients",index_col=0)
def import_biosamples():
    return pd.read_excel(SUPPLEMENTARY_TABLES_PATH,sheet_name="2. Biosamples",index_col=0)
def import_amplicons():
    return pd.read_excel(SUPPLEMENTARY_TABLES_PATH,sheet_name="3. Amplicons")
def import_genes():
    return pd.read_excel(SUPPLEMENTARY_TABLES_PATH,sheet_name="5. Gene amplifications",
                         na_values = ['unknown'],
                         converters={'gene_cn': float, 'is_canonical_oncogene': bool})