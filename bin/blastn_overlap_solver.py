#!/usr/bin/env python3

#############################################################
# Created by Jorge R Grande - HUMV - Santander
# Modified by Rosalía Palomino-Cabrera
# blastn_overlap_solver.py takes an .tsv table from blastn
# and outputs a parsed table without overlapping hits.
# Best hit is chosen based on Blast score

import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def get_arguments():

    parser = argparse.ArgumentParser(prog = 'blastn_overlap_solver.py', description = 'Solves overlapping hits from Blastn search.')

    input_group = parser.add_argument_group('Input')
    input_group.add_argument('-i', '--input_file', dest="input_file", required=True, help="Required. Input file to be filtered", type=os.path.abspath)

    output_group = parser.add_argument_group('Output')
    output_group.add_argument('-o', '--out_dir', dest='out_dir', required=True, help='Required. Final output folder for reports', type=os.path.abspath)

    arguments = parser.parse_args()

    return arguments

def main():
    args = get_arguments()

    # Parameters
    report = args.input_file
    out_dir = args.out_dir
    outfile = Path(report).stem + '_filtered.tsv'
    report_out = os.path.join(out_dir, outfile)

    # Input processing
    df = pd.read_table(report, sep='\t')
    header = ['IS', 'contig', 'qstart', 'qend', 'qlen', 'sstart0', 'send0', 'slen', 'pident', 'qcovhsp', 'length', 'mismatch', 'score', 'evalue']
    df.columns = header
    df['contig'] = df['contig'].str.replace(r'^[a-zA-Z]+_', '', regex=True).astype(int)

    coords = np.sort(df[['sstart0', 'send0']], axis=1)  # Redefines sequence coordinates in negative frames
    df['sstart'] = coords[:, 0]
    df['send'] = coords[:, 1]
    df.drop(columns=['sstart0', 'send0'], inplace=True)

    # Drops hits with identical coords, keeps higher score hit
    df.sort_values(by=["contig", "send","score"], ascending=[True, True, False], inplace=True)
    df.drop_duplicates(subset=["contig", 'send'], keep='first', inplace=True)           # End coords
    df.sort_values(by=["contig", "sstart","score"], ascending=[True, True, False], inplace=True)
    df.drop_duplicates(subset=["contig", 'sstart'], keep='first', inplace=True)         # Start coords
    df.reset_index(drop=True, inplace=True)

    # Drops overlapping hits
    filtered = []
    for _, row in df.iterrows():
        # Checks hits are in the same contig and the current IS starts before the previous one ends
        while filtered and row["contig"] == filtered[-1]['contig'] and row['sstart'] < filtered[-1]['send']:
            if row['score'] > filtered[-1]['score']:
                filtered.pop()
            else:
                break
        else:
            filtered.append(row)
    filtered_df = pd.DataFrame(filtered)

    # Output
    filtered_df.to_csv(report_out, sep='\t', index=False)

if __name__ == '__main__':
    main()
