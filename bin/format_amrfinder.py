#!/usr/bin/env python3

"""
Script to reformat AMRFINDERPLUS output, mimicking Abricate .txt outfile.
"""

import argparse
import pandas as pd
import numpy as np

def load_tab(input):
    df = pd.read_table(input, header=0)
    df['Start'] = df['Start'].astype(int)
    df['Stop'] = df['Stop'].astype(int)
    return df[['Contig id', 'Start', 'Stop',
               'Strand', 'Element symbol',
               '% Coverage of reference',
               '% Identity to reference', 'Scope',
               'Closest reference accession',
               'Element name', 'Subclass']]

def reformat_tab(df):
    reformed_df = df.rename(columns={'Contig id': 'SEQUENCE',
                                     'Start': 'START',
                                     'Stop': 'END',
                                     'Strand': 'STRAND',
                                     'Element symbol': 'GENE',
                                     '% Coverage of reference': '%COVERAGE',
                                     '% Identity to reference': '%IDENTITY',
                                     'Scope': 'DATABASE',
                                     'Closest reference accession': 'ACCESSION',
                                     'Element name': 'PRODUCT',
                                     'Subclass': 'RESISTANCE'}
                                     )
    return reformed_df

def main(input, output):
    df_amr = load_tab(input)
    reformed_df = reformat_tab(df_amr)
    reformed_df.to_csv(output, index=False, sep='\t')
    print(f"Wrote reformated output to {output}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reformat AMRFinderPlus.tsv file to mimick Abricate.txt output.')
    parser.add_argument('-i', '--input', required=True, help='Path to input AMRFinderPlus file')
    parser.add_argument('-o', '--output', required=True, help='Path to output reformated file')
    args = parser.parse_args()

    main(args.input, args.output)
