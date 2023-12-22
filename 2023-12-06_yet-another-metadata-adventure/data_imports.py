#! python3
# Library of functions to load metadata from various places, and merge into unified patient, biosample tables.
# Usage: generate_biosample_table()
# Requires: 
# 
# Owen Chapman

import pandas as pd
import pathlib
import os

## Function to load metadata from the AmpliconClassifier results
## Get this file from /expanse/lustre/projects/csd677/collab/projects/pedpancan/AmpliconClassifier/batch/inputs
def get_pedpancan_biosamples_from_AC(include_x01=False,path='../data/source/AmpliconClassifier/pedpancan_summary_map.txt'):
    path = pathlib.Path(path)
    df = pd.read_csv(path, sep='\t', header=None, index_col=0, names = ["biosample","file"])
    if not include_x01:
        df['firstletter']=df.index.map(lambda x: x[0])
        df = df[df.firstletter != 'P']
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
    # remove suffix from pnoc sample ids
    df.sample_id = df.sample_id.map(lambda x: '-'.join(x.split('-')[:2]) if x.startswith("7316-") else x)
    
    df = df.replace({
        'Tumor Descriptor':{
            "initial CNS Tumor": "Diagnosis",
            "Not Applicable":None,
            "Unavailable":None,
            "Initial CNS Tumor": "Diagnosis",
            "Progressive Disease Post-Mortem":"Progressive",
        },
        'gender':{
            "Not Reported":None
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
def import_cbtn_biosample_metadata(include_X01=False):
    if include_X01:
        df = pd.concat([import_x00_biosample_metadata(),import_x01_biosample_metadata(),import_pnoc_biosample_metadata()])
    else:
        df = pd.concat([import_x00_biosample_metadata(),import_pnoc_biosample_metadata()])
    cohort = get_pedpancan_biosamples_from_AC()
    df = df[df.index.isin(cohort)]
    df = clean_cavatica_biosample_metadata(df)
    return df

## Functions to open & preprocess opentarget histology data.
## Get histologies.tsv from https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/analyses/molecular-subtyping-integrate/results/histologies.tsv
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
def import_opentarget_histologies_files(path='../data/local/opentarget/histologies.tsv'):
    '''
    Get this file from /Users/ochapman/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/projects/2023-pedpancan/data/opentarget/histologies.tsv
    '''
    path = pathlib.Path(path)
    df = pd.read_csv(path,sep='\t',index_col=0)
    df = clean_opentarget_histologies_files(df)
    return df
def import_pedcbioportal_metadata(path="../data/source/pedcbioportal/openpbta-biosample-metadata.tsv"):
    path = pathlib.Path(path)
    df = pd.read_csv(path, sep='\t',index_col=0)
    return df
def get_cbtn_cell_lines():
    df = import_pedcbioportal_metadata()
    df = df[df.SAMPLE_TYPE == "Derived Cell Line"]
    return df.SPECIMEN_ID.str.cat(sep=';').split(';')

#duplicated_cbtn_samples = ["BS_DRVEFVQ5","BS_169P1QCA","BS_QG6V29H7",
#                           "BS_6GV08HTE","BS_B4PPG3X5","BS_2X60Q1ET","BS_S791VC80","BS_STNH7YSX","BS_3Z40EZHD","BS_ZR75EKKX",
#                           "BS_KQPCYZ2K","BS_Z64NEPNE","BS_KQRAHH6Y","BS_KH3859M5","BS_HJ7HYZ7N","BS_AH3RVK53",
#                           "BS_5S8VXASX","BS_JEZBA2EW","BS_XNYQS1WG","BS_P99S85CY",
#                           "BS_M5FM63EB","BS_M0B42FPR","BS_TX2WGF8K","BS_RENPFFNK","BS_R6CKWZW6",
#                           "BS_ZSH09N84","BS_6TMPZKSZ","BS_B91XGSA5","BS_XQF18WZP","BS_0TCRV9AC",
#                           "BS_2J4FG4HV",
#                           "BS_EE73VE7V","BS_5968GBGT","BS_BQ81D2BP","BS_3VKW5988", # duplicate samples from PT_KTRJ8TFY autopsy
#                           "BS_AK9BV52G","BS_X5VN0FW0","BS_D6STCMQS","BS_22VCR7DF","BS_1Q524P3B" # duplicate samples from PT_KZ56XHJT autopsy
#                          ]
#nontumor_samples = ["BS_MCM78YPC","BS_886M7JMG","BS_TPX7YY57"] # Epilepsy, Arteriovenous malformation, and Reactive connective tissue respectively

def propagate(df,dest,source,rename=False):
    '''
    Replace NA values in dest with those in source, then drop source and rename dest
    '''
    df[dest].fillna(df[source], inplace=True)
    df.drop(source, inplace=True, axis=1)
    if rename:
        df = df.rename(columns={dest:rename})
    return df
def consensus(df,dest,source,rename=False):
    '''
    Check that values in dest and source agree, if not then set to NA. Drop source and rename dest.
    '''
    df.loc[df[dest] != df[source], dest] = None
    df.drop(source, inplace=True, axis=1)
    if rename:
        df = df.rename(columns={dest:rename})
    return df    

## Integrate all CBTN data available
def generate_cbtn_biosample_table(verbose=0):
    '''
    Generate a metadata table of cbtn biosamples.
    verbose:
        0: most useful metadata only. These are included in the Supplmentary Table.
        1: includes some extra columns. Useful for generating the patient metadata table.
        2: includes lots of extra columns
    '''
    df = pd.DataFrame(index=get_pedpancan_biosamples_from_AC())
    cavatica_data = import_cbtn_biosample_metadata()
    df = pd.merge(left=df,how='inner',right=cavatica_data,left_index=True,right_index=True)
    opentarget_data = import_opentarget_histologies_files()
    df = pd.merge(left=df,how='left',right=opentarget_data,left_index=True,right_index=True,suffixes=(None,"_y"))

    # For CAVATICA annotations which are missing, propagate those from opentarget.
    df = propagate(df,"primary_site","primary_site_y")
    df = propagate(df,"age_at_diagnosis","age_at_diagnosis_days")
    df = propagate(df,"Tumor Descriptor","tumor_descriptor")
    df = consensus(df,"gender","reported_gender","sex")
    
    # Rename columns
    df.index.name = "biosample_id"
    df = df.rename(columns={
        'Kids First Participant ID':'patient_id',
        'Tumor Descriptor':'tumor_history',
        'case_id':'external_patient_id',
        'sample_id':'external_sample_id',
    })

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
    # Drop nontumor samples
    #df = df.drop(nontumor_samples)
    
    # Drop cell lines
    #df = df[~df.index.isin(get_cbtn_cell_lines())]
    
    # Mark duplicates
    #df['in_deduplicated_sample_cohort'] = True
    #df.loc[duplicated_cbtn_samples,'in_deduplicated_sample_cohort'] = False
    return df

## SJ data
# duplicated_sj_samples are biosamples where we suspect that the same tumor was sequenced twice at the same timepoint. 
# These are flagged and removed arbitrarily.
#duplicated_sj_samples = ['SJST030043_D1','SJST030131_D1','SJMEL001003_D2','SJOS001115_D1','SJMB009_E',
#                         'SJDSRCT030041_D3','SJBT030081_D2','SJWLM030180_D2','SJHGG017_D','SJEWS030228_D2',
#                         'SJST030383_D1','SJLGG017_D','SJLGG030611_D2','SJOS030876_D2','SJOS001101_M1','SJHM030702_D2']
# what to do with SJOS001101_M1?

def import_sj_sample_info(path="../data/local/sjcloud/SAMPLE_INFO_2022-03-02.tsv"):
    path = pathlib.Path(path)
    df = pd.read_csv(path,sep='\t',index_col="sample_name")
    return df
def clean_sj_biosample_metadata(df):
    '''
    Clean known errors in the sj metadata, and unify ontologies, units etc.
    '''
    df = df.replace({
        'attr_age_at_diagnosis':{
            "Not Available": None
        },
        'attr_sex':{
            "Not Available":None
        }
    })
    
    # Convert age from years to days
    df['attr_age_at_diagnosis'] = (pd.to_numeric(df['attr_age_at_diagnosis'],errors='coerce')*365.25).round()
    
    return df
def generate_sj_biosample_table(verbose=0):
    '''
    Notes:
    sj_diseases != attr_oncotree_disease_code = sj_associated_diagnoses_disease_code
    attr_diagnosis != sj_long_disease_name != sj_associated_diagnoses
    '''
    df = pd.DataFrame(index=get_pedpancan_biosamples_from_AC())
    columns = ['subject_name','sample_type','attr_age_at_diagnosis','attr_sex','sj_long_disease_name','sj_diseases','attr_oncotree_disease_code','sj_dataset_accessions']
    add = import_sj_sample_info()
    add = clean_sj_biosample_metadata(add)
    add = add[(add.sequencing_type == 'WGS') & (add.file_type == 'BAM') & add.file_path.str.endswith('.bam')]
#    add = add.sort_values(columns)
#    add = add.loc[~add.index.duplicated()]
    df = pd.merge(left=df,how='inner',right=add, left_index=True, right_index=True)
    
    # Rename columns
    df.index.name = "biosample_id"
    df = df.rename(columns={
        'subject_name':'patient_id',
        'sample_type':'tumor_history',
        'attr_sex':'sex',
        'sj_dataset_accessions':'cohort',
        'attr_age_at_diagnosis':'age_at_diagnosis',
    })
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
    
    # Mark duplicates
    #df['in_deduplicated_sample_cohort'] = True
    #df.loc[duplicated_sj_samples,'in_deduplicated_sample_cohort'] = False
    return df

## Generate the unified biosample table
def unify_tumor_diagnoses(df):
    '''
    TODO: Add code mapping SJ and CBTN tumor types to Sunita's.
    '''
    # @Rishaan Add classification code here
    
    # drop tumor type annotations now that we have a unified diagnosis.
    df = df.drop(["disease_type","dkfz_v11_methylation_subclass","dkfz_v11_methylation_subclass_score",
                  "dkfz_v12_methylation_subclass","dkfz_v12_methylation_subclass_score","molecular_subtype","harmonized_diagnosis",
                  "broad_histology","short_histology","sj_long_disease_name","sj_diseases"
                 ],axis=1)
    return df






## Annotate with ecDNA status
def annotate_with_ecDNA(df,path="../data/Supplementary Tables.xlsx"):
    '''
    Annotate biosamples with ecDNA status.
    Inputs:
        df: pd.DataFrame. Must be indexed by biosample.
        path: path to AmpliconClassifier results.
    '''
    # load AC results
    if path.endswith("Supplementary Tables.xlsx"):
        ac = pd.read_excel(path,sheet_name="3. Amplicons")
    else:
        ac = pd.read_excel(path,index_col=0)
    
    # Aggregate by biosample
    ac_agg = ac.groupby("sample_name").sum().ecDNA_amplicons
    df = df.join(ac_agg)
    df = df.rename(columns={"ecDNA_amplicons":"ecDNA_sequences_detected"})
    df["ecDNA_sequences_detected"].fillna(0,inplace=True)
    return df

def generate_biosample_table():
    df = pd.concat([generate_cbtn_biosample_table(),generate_sj_biosample_table()])
    df = unify_tumor_diagnoses(df)
    df = annotate_with_ecDNA(df)
    return df

## Imports for Sunita's data
def import_sunita_classifications(path='../data/combinedamplicons.xlsx'):
    df = pd.read_excel(path)
    df = df[["sample_ID","cancer_type"]]
    df = df.drop_duplicates()
    df = df.set_index("sample_ID")
    return df
