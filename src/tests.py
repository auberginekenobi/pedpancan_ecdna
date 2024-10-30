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
        test_dubois_patient_integration(p),
        test_dubois_subtype_integration(b),
        test_patient_amp_class(p,b),
        test_biosample_amp_class(b,a),
    ])
    for r in results:
        print(r)

    print("passed all tests!")
    return

