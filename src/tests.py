import inspect
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
    assert ibset == ipset, text_venn2(ibset, ipset)

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
    assert iaset == ibset, text_venn2(naset, nbset)
    
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
        assert(set_d | set_w == set())
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

def run_all_tests(patients = None, biosamples = None, amplicons = None):
    # Generate tables once
    if biosamples is None:
        biosamples = generate_biosample_table()
    if patients is None:
        patients = generate_patient_table(biosamples)
    if amplicons is None:
        amplicons = generate_amplicon_table(biosamples)
    p,b,a = patients,biosamples,amplicons

    # Run tests
    results = (r for r in [
        test_patient_amp_class(p,b),
        test_biosample_amp_class(b,a),
        test_dubois_patient_integration(p),
        test_dubois_subtype_integration(b),
        test_dubois_cancer_type_disambiguation(b),
        test_sample_deduplication_max_ecDNA(b)
    ])
    for r in results:
        print(r)

    print("passed all tests!")
    return

