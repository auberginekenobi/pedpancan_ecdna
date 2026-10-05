import inspect
import numpy as np
import pandas as pd
import warnings
from data_imports import *


def text_venn2(s1, s2):
    print(f'size of set 1: {len(s1)}')
    print(f'size of set 2: {len(s2)}')
    print(f'Samples in s1 not in s2: {len(s1 - s2)}')
    print(f'Samples in s2 not in s1: {len(s2 - s1)}')
    print(f'Overlap: {len(s1 & s2)}')

def test_patient_amp_class(patients=None,biosamples=None):
    if patients is None:
        patients = generate_patient_table()
    if biosamples is None:
        biosamples = generate_biosample_table()

    ebset = set(biosamples[biosamples.amplicon_class == 'ecDNA'].patient_id)
    epset = set(patients[patients.amplicon_class == 'ecDNA'].index)
    assert ebset == epset, text_venn2(ebset, epset)

    ibset = set(biosamples[biosamples.amplicon_class == 'intrachromosomal'].patient_id) - ebset
    ipset = set(patients[patients.amplicon_class == 'intrachromosomal'].index)
    assert ibset == ipset, text_venn2(ibset, ipset)

    nbset = set(biosamples[biosamples.amplicon_class == 'no amplification'].patient_id) - ebset - ibset
    npset = set(patients[patients.amplicon_class == 'no amplification'].index)
    assert nbset == npset, text_venn2(nbset, npset)

    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_biosample_amp_class(biosamples=None,amplicons=None):
    if biosamples is None:
        biosamples = generate_biosample_table()
    if amplicons is None:
        amplicons = generate_amplicon_table()
    
    ebset = set(biosamples[biosamples.amplicon_class == 'ecDNA'].index)
    easet = set(amplicons[amplicons['ecDNA+'] == 'Positive']['sample_name'])
    assert easet == ebset, text_venn2(easet, ebset)
    
    ibset = set(biosamples[biosamples.amplicon_class == 'intrachromosomal'].index)
    iaset = set(amplicons[(amplicons['BFB+'] == 'Positive') |
                           (amplicons['amplicon_decomposition_class'].isin(['Complex-non-cyclic','Linear']))
                ]['sample_name']) - easet
    assert iaset == ibset, text_venn2(iaset, ibset)

    nbset = set(biosamples[biosamples.amplicon_class == 'no amplification'].index)
    naset = set(amplicons.sample_name) - iaset - easet
    assert nbset >= naset, text_venn2(nbset, naset)
    
    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_dubois_subtype_integration(biosamples = None):
    if biosamples is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            biosamples = generate_biosample_table()
        assert(len(biosamples[biosamples.cancer_subclass == 'HGG_H3K27']) > 0)
    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_dubois_patient_integration(patients = None):
    if patients is None:
        patients = generate_patient_table()
    val = patients.loc['SJ000101','OS_status']
    assert val == 'Deceased', val # was NA without Dubois data
    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_dubois_cancer_type_disambiguation(biosamples = None):
    '''
    In SJ annotations, WT means Wilms' tumor. In Dubois, it means H3 wild-type. 
    Added code to disambiguate based on source ontology.
    '''
    if biosamples is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            biosamples = generate_biosample_table()
        assert(biosamples.loc['SJHGG017_D','cancer_type'] == 'HGG')
        set_d = set(import_dubois_supplementary_data().index)
        set_w = set(biosamples[biosamples.cancer_type == 'WLM'].index)
        assert(set_d & set_w == set())
    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_sample_deduplication_max_ecDNA(biosamples = None):
    '''
    Sample deduplications should take the sample with the most ecDNA amps to to ameliorate the 
    problem of ecDNAs missing in downstream analyses.
    Eg. for PT_XA98HG1C, BS_5JC116NM has 1 ecDNA but BS_W37QBA12 and BS_2J4FG4HV have 2,
    so BS_W37QBA12 or BS_2J4FG4HV should be the deduplicated sample.
    '''
    if biosamples is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            biosamples = generate_biosample_table()
    assert not biosamples.loc['BS_5JC116NM','in_unique_tumor_set']
    assert not biosamples.loc['BS_5JC116NM','in_unique_patient_set']
    assert biosamples.loc['BS_W37QBA12','in_unique_tumor_set'] or biosamples.loc['BS_2J4FG4HV','in_unique_tumor_set']
    assert biosamples.loc['BS_W37QBA12','in_unique_patient_set'] or biosamples.loc['BS_2J4FG4HV','in_unique_patient_set']
    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_all_cancer_types_annotated(biosamples=None):
    if biosamples is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            biosamples = generate_biosample_table()
    assert sum(biosamples.cancer_type.isna()) == 0
    return f'pass: {inspect.currentframe().f_code.co_name}'

def read_consent_withdrawals(file='../../data/source/sjcloud/RTCG_sample_removal_202609_ic.xlsx'):
    return pd.read_excel(file)

def test_consent_withdrawn(patients=None,biosamples=None):
    if biosamples is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            biosamples = generate_biosample_table()
    if patients is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            patients = generate_patient_table()
    withdrawals = read_consent_withdrawals()
    pt_mask = patients.index.isin(withdrawals.subject_name)
    try:
        assert not bool(pt_mask.any())
    except AssertionError:
        print(f'Patients should be excluded: {patients.index[pt_mask].tolist()}'); raise
    bs_mask = biosamples.index.isin(withdrawals.sample_name)
    try:
        assert not bool(bs_mask.any())
    except AssertionError:
        print(f'Biosamples should be excluded: {biosamples.index[bs_mask].tolist()}'); raise
    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_drop_cell_lines(biosamples = None):
    def get_cell_lines_from_opentarget(path='../../data/source/opentarget/histologies.tsv',verbose=False):
        path = pathlib.Path(path)
        df = pd.read_csv(path,sep='\t',index_col=0,low_memory=False)
        df = df[(df.composition == 'Derived Cell Line') & (df.index.str.startswith('BS'))]
        return df.index.unique().tolist()
    cell_lines = get_cell_lines_from_opentarget()
    if biosamples is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            biosamples = generate_biosample_table()
    bs_mask = biosamples.index.isin(cell_lines)
    try:
        assert not bool(bs_mask.any())
    except AssertionError:
        print(f'Cell lines should be excluded: {biosamples.index[bs_mask].tolist()}'); raise
    return f'pass: {inspect.currentframe().f_code.co_name}'
    

def test_drop_misc_patients(patients = None):
    drop_ids = [
        'PT_AQ2Q3JMC', # patient assigned multiple PT_ids
        'PT_EDG0Q7P4' # patient assigned multiple PT_ids
    ]
    if patients is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            patients = generate_patient_table()
    pt_mask = patients.index.isin(drop_ids)
    try:
        assert not bool(pt_mask.any())
    except AssertionError:
        print(f'Patients should be excluded: {patients.index[pt_mask].tolist()}'); raise
    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_drop_misc_biosamples(biosamples = None):
    drop_ids = [
        'BS_HJ7HYZ7N', # mis-annotated normal
    ]
    if biosamples is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            biosamples = generate_biosample_table()
    bs_mask = biosamples.index.isin(drop_ids)
    try:
        assert not bool(bs_mask.any())
    except AssertionError:
        print(f'Biosamples should be excluded: {biosamples.index[bs_mask].tolist()}'); raise
    return f'pass: {inspect.currentframe().f_code.co_name}'

def assert_in_range(df, column, low=-np.inf, high=np.inf, inclusive='both'):
    '''
    Assert that all non-missing values of df[column] lie between low and high; print offending rows otherwise.
    '''
    values = pd.to_numeric(df[column])
    mask = values.notna() & ~values.between(low, high, inclusive=inclusive)
    try:
        assert not bool(mask.any())
    except AssertionError:
        print(f'{column} should be between {low} and {high} (inclusive={inclusive}): {values[mask].to_dict()}'); raise

def test_age_range(patients = None, biosamples = None):
    if biosamples is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            biosamples = generate_biosample_table()
    if patients is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            patients = generate_patient_table()
    assert_in_range(biosamples, 'age_at_diagnosis', 0, 36525)
    assert_in_range(biosamples, 'age_at_collection', 0, 36525)
    assert_in_range(patients, 'age_at_diagnosis', 0, 36525)
    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_tumor_fraction_range(patients = None, biosamples = None):
    if biosamples is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            biosamples = generate_biosample_table()
    if patients is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            patients = generate_patient_table()
    assert_in_range(biosamples, 'tumor_fraction_THetA2', 0, 1, inclusive='neither')
    assert_in_range(patients, 'tumor_fraction_THetA2', 0, 1, inclusive='neither')
    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_survival_range(patients = None):
    if patients is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=UserWarning)
            patients = generate_patient_table()
    assert_in_range(patients, 'OS_months', high=365.25*12, inclusive='neither')
    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_ecDNA_amplicons_range(amplicons = None):
    if amplicons is None:
        amplicons = generate_amplicon_table()
    assert_in_range(amplicons, 'ecDNA_amplicons', low=0)
    return f'pass: {inspect.currentframe().f_code.co_name}'

def test_gene_cn_range(genes = None):
    if genes is None:
        genes = generate_gene_table()
    assert_in_range(genes.replace({'gene_cn': {'unknown': np.nan}}), 'gene_cn', 0, 1000) # AmpliconClassifier reports some gene_cn as 'unknown'
    return f'pass: {inspect.currentframe().f_code.co_name}'

def run_all_tests(patients = None, biosamples = None, amplicons = None, genes = None):
    # Generate tables once
    if biosamples is None:
        biosamples = generate_biosample_table()
    if patients is None:
        patients = generate_patient_table(biosamples)
    if amplicons is None:
        amplicons = generate_amplicon_table(biosamples)
    if genes is None:
        genes = generate_gene_table(biosamples)
    p,b,a,g = patients,biosamples,amplicons,genes

    # Run tests
    results = (r for r in [
        test_patient_amp_class(p,b),
        test_biosample_amp_class(b,a),
        test_dubois_patient_integration(p),
        test_dubois_subtype_integration(b),
        test_dubois_cancer_type_disambiguation(b),
        test_sample_deduplication_max_ecDNA(b),
        test_all_cancer_types_annotated(b),
        test_consent_withdrawn(p,b),
        test_drop_cell_lines(b),
        test_drop_misc_biosamples(b),
        test_drop_misc_patients(p),
        test_age_range(p,b),
        test_tumor_fraction_range(p,b),
        test_survival_range(p),
        test_ecDNA_amplicons_range(a),
        test_gene_cn_range(g)
    ])
    for r in results:
        print(r)

    print("passed all tests!")
    return

