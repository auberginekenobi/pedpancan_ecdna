#! python3
# Library of functions to assemble RNA-seq library metadata for the gene expression analyses.
#
# Owen Chapman

import pathlib
import pandas as pd
import data_imports

def resolve_sj_duplicate_libraries(df,fc_dir='../../data/source/gene-expression/sjcloud/FEATURE_COUNTS'):
    '''
    Some SJ biosamples have more than one RNA libraries with the same name. Disambiguate by file_size.
    Inputs:
        df: pd.DataFrame with one row per file_id, and rna_id and file_size columns.
        fc_dir: directory of <biosample>.RNA-Seq.feature-counts.txt files.
    '''
    dup = df.rna_id.duplicated(keep=False)
    held_size = df.rna_id[dup].map(lambda x: pathlib.Path(fc_dir,f'{x}.RNA-Seq.feature-counts.txt').stat().st_size).reindex(df.index)
    resolved = df[~dup | (df.file_size == held_size)]
    unresolved = set(df.rna_id) - set(resolved.rna_id)
    if unresolved or resolved.rna_id.duplicated().any():
        raise ValueError(f'Could not identify exactly one held file for SJ libraries: {sorted(unresolved | set(resolved.rna_id[resolved.rna_id.duplicated()]))}')
    return resolved

def import_sj_rna_metadata(biosamples_tbl=None,path='../../data/source/sjcloud/SAMPLE_INFO_20260901.txt'):
    '''
    Return metadata for SJ RNA-seq libraries in our expression matrix from S.T.2 biosamples, indexed by rna_id (matrix column).
    SJ RNA-seq is indexed by event, so each library's event is its own biosample.
    Inputs:
        biosamples_tbl: S.T.2. If None, import it.
        path: SJ Cloud SAMPLE_INFO listing the RNA-seq feature counts files.
    '''
    if biosamples_tbl is None:
        biosamples_tbl = data_imports.import_biosamples()
    held = data_imports.get_sj_rna_biosample_ids()
    df = pd.read_csv(pathlib.Path(path),sep='\t',low_memory=False)
    df = df[(df.sequencing_type == 'RNA-Seq') & (df.file_type == 'FEATURE_COUNTS')]
    df = df[df.sample_name.isin(held) & df.sample_name.isin(biosamples_tbl.index)]
    df = df[df.attr_tissue_preservative != 'FFPE'] # FFPE degrades RNA; keep 'Fresh/Frozen(NGS), FFPE(other)' values.
    df = df.rename(columns={
        'sample_name':'rna_id',
        'attr_library_selection_protocol':'rna_selection',
        'attr_lab_strandedness':'strandedness',
        'attr_read_length':'read_length',
    })
    df = resolve_sj_duplicate_libraries(df)
    df['cohort'] = 'SJC'
    df['patient_id'] = df.rna_id.map(biosamples_tbl.patient_id) # SJ subject_name has changed between releases for 9 SJACT cases
    df['event_id'] = df.rna_id
    df = df[['rna_id','cohort','patient_id','event_id','rna_selection','strandedness','read_length']]
    df = df.replace({
        'rna_selection':{
            'Not Available':pd.NA,
        },
        'strandedness':{
            'Stranded-Reverse':'Stranded',
            'Not Available':pd.NA,
        }
    })
    # Total RNA is always reverse stranded and PolyA always unstranded; infer the missing annotation from the other.
    df.loc[(df.rna_selection == 'Total').fillna(False) & df.strandedness.isna(),'strandedness'] = 'Stranded'
    df.loc[(df.strandedness == 'Stranded').fillna(False) & df.rna_selection.isna(),'rna_selection'] = 'Total'
    df = df[df.rna_selection.notna()] # drop library prep unknown
    df['read_length'] = df.read_length.map({'100':'101','101':'101','125':'126','126':'126'}).fillna('NEC').astype('string') # read length per read; other values are singletons, ambiguous or uninterpretable
    df = df.set_index('rna_id')
    return df

def import_cbtn_rna_metadata(biosamples_tbl=None,path='../../data/source/opentarget/histologies_v15.tsv'):
    '''
    Return metadata for CBTN RNA-seq libraries in our expression matrix from S.T.2 collection events, indexed by rna_id (matrix column).
    RNA-seq has its own biospecimen and aliquot IDs, so libraries are linked to S.T.2 biosamples by collection event (7316- sample id).
    Inputs:
        biosamples_tbl: S.T.2. If None, import it.
        path: OpenPedCan histologies matching the expression matrix release.
    '''
    if biosamples_tbl is None:
        biosamples_tbl = data_imports.import_biosamples()
    held = data_imports.get_cbtn_rna_biospecimen_ids()
    df = pd.read_csv(pathlib.Path(path),sep='\t',low_memory=False)
    df = df[df.Kids_First_Biospecimen_ID.isin(held)]
    df = df[(df.cohort == 'PBTA') & (df.experimental_strategy == 'RNA-Seq') & (df.sample_type == 'Tumor') & (df.composition == 'Solid Tissue')]
    df = df[df.RNA_library.isin(['poly-A','poly-A stranded','stranded'])] # drop exome capture
    drop_center = ["The Translational Genomics Research Institute","Broad Institute"]
    df = df[~df.seq_center.isin(drop_center)]
    df['event_id'] = df.sample_id.map(lambda x: '-'.join(x.split('-')[:2]) if x.startswith("7316-") else x) # remove suffix from pnoc sample ids, as in S.T.2
    event_patient = biosamples_tbl.dropna(subset=['external_sample_id']).groupby('external_sample_id').patient_id.first()
    df = df[df.event_id.isin(event_patient.index)]
    df = df.rename(columns={
        'Kids_First_Biospecimen_ID':'rna_id',
    })
    df['cohort'] = 'CBTN'
    df['patient_id'] = df.event_id.map(event_patient)
    df['rna_selection'] = df.RNA_library.map({'poly-A':'PolyA','poly-A stranded':'PolyA','stranded':'Total'})
    df['strandedness'] = df.RNA_library.map({'poly-A':'Unstranded','poly-A stranded':'Stranded','stranded':'Stranded'})
    df = df.replace({
        'seq_center':{
            'BGI@CHOP Genome Center':'BGI@CHOP',
            'BGI':'BGI@CHOP',
            'Not Available':pd.NA
        },
    })
    
    df = df[['rna_id','cohort','patient_id','event_id','rna_selection','strandedness','seq_center']]
    df = df.set_index('rna_id')
    return df

def generate_rna_metadata_table(biosamples_tbl=None):
    '''
    Return metadata for the RNA-seq libraries we analyze, one row per library, indexed by rna_id (expression matrix column).
    Event-level annotations from S.T.2 are joined by event_id. Where WGS biosamples from the same event disagree, take the value
    from the event's representative biosample (in_unique_tumor_set).
    Inputs:
        biosamples_tbl: S.T.2. If None, import it.
    '''
    if biosamples_tbl is None:
        biosamples_tbl = data_imports.import_biosamples()
    df = pd.concat([import_sj_rna_metadata(biosamples_tbl),import_cbtn_rna_metadata(biosamples_tbl)])
    if not df.index.is_unique:
        raise ValueError(f'Duplicate rna_ids: {df.index[df.index.duplicated()].tolist()}')
    events = biosamples_tbl.assign(event_id=biosamples_tbl.external_sample_id.fillna(pd.Series(biosamples_tbl.index,index=biosamples_tbl.index)))
    events = events.sort_values('in_unique_tumor_set',ascending=False,kind='stable') # representative biosample first
    # TODO: event 7316-3217 has WGS biosamples annotated as both Diagnosis (BS_H8NWA41N) and Progressive (BS_0ATJ22QA).
    events = events.groupby('event_id').agg(
        cancer_type=('cancer_type','first'),
        cancer_subclass=('cancer_subclass','first'),
        tumor_history=('tumor_history','first'),
        sex=('sex','first'),
        age_at_collection=('age_at_collection','first'),
        age_at_diagnosis=('age_at_diagnosis','first'),
        ecDNA_any=('ecDNA_sequences_detected',lambda x: (x > 0).any()),
    )
    df = df.join(events,on='event_id')
    return df