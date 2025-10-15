#!/usr/bin/env python3

# Write a set of cancer genes to a file
# Input: path of the COSMIC Cancer Gene Census file, eg. Cosmic_CancerGeneCensus_v101_GRCh38.tsv.
# Output: path to your new gene list, eg. cancer_genes.txt.

import argparse
import pandas as pd

def import_cgc(file):
    """Read COSMIC Cancer Gene Census TSV file."""
    return pd.read_table(file)

def get_cancer_genes(file, drop_fusions=True):
    """Extract cancer gene symbols from the CGC file."""
    df = import_cgc(file)
    df = df.dropna(subset=['ROLE_IN_CANCER'])
    if drop_fusions:
        df = df[df['ROLE_IN_CANCER'] != 'fusion']
    return set(df['GENE_SYMBOL'])

def main():
    parser = argparse.ArgumentParser(
        description="Extract gene symbols from COSMIC Cancer Gene Census TSV file."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to COSMIC Cancer Gene Census TSV file"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output path for text file containing one gene symbol per line"
    )

    args = parser.parse_args()

    genes = get_cancer_genes(args.input)

    with open(args.output, "w") as f:
        for gene in sorted(genes):
            f.write(f"{gene}\n")

if __name__ == "__main__":
    main()
