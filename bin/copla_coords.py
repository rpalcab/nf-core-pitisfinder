#!/usr/bin/env python3

"""
Script to retrieve plasmid biomarkers from COPLA output.
"""

import argparse
import pandas as pd
import os


def read_mob_file(mob_path):
    if os.path.getsize(mob_path) == 0:
        return pd.DataFrame(columns=['#Hit_Id', 'Gene', 'Sequence_coverage', 'i-evalue', 'Reference_system'])

    df_raw = pd.read_csv(mob_path, sep="\t", header=None)
    df = df_raw.iloc[:, [0, 1, 3, 6]].copy()
    df.columns = ['#Hit_Id', 'Gene', 'Sequence_coverage', 'i-evalue']
    df['Reference_system'] = "MOB"
    return df


def read_rep_file(rep_path):
    df = pd.read_csv(rep_path, sep="\t", header=0, usecols=[1, 2, 4, 5])
    df.rename(columns={'Plasmid': 'Gene'}, inplace=True)

    if not df.empty:
        df[['start', 'end']] = df['Position in contig'].str.split(r'\.\.', expand=True)
        df.drop(columns=['Position in contig'], inplace=True)
    else:
        df[['start', 'end']] = "", ""

    df['strand'] = pd.NA
    df['Coverage'] = pd.NA
    df['i-evalue'] = pd.NA
    df['Reference_system'] = "Replicon"
    return df


def read_conj_file(conj_path):
    if conj_path and os.path.getsize(conj_path) > 0:
        df = pd.read_csv(conj_path, sep="\t", header=0, usecols=[0, 4, 5, 10, 13])
        df.columns = ['#Hit_Id', 'Gene', 'Reference_system', 'i-evalue', 'Sequence_coverage']

        # Replace values conditionally
        df.loc[df['Reference_system'].str.startswith("type", na=False), 'Reference_system'] = 'MPF'
        df.loc[df['Reference_system'] == 'CONJ', 'Reference_system'] = 'MOB'
        return df
    return pd.DataFrame(columns=['#Hit_Id', 'Gene', 'Reference_system', 'i-evalue', 'Sequence_coverage'])


def main(gene_pos_path, mob_path, rep_path, conj_path, output_path):
    # Read gene positions
    df_pos = pd.read_csv(gene_pos_path, sep="\t", header=0)

    # Read MOB, CONJ, REP
    df_mob = read_mob_file(mob_path)
    df_conj = read_conj_file(conj_path)
    df_rep = read_rep_file(rep_path)

    # Combine MOB and CONJ
    df_all = pd.concat([df_mob, df_conj], ignore_index=True)

    # Merge with gene positions
    df_all = df_all.merge(df_pos, left_on='#Hit_Id', right_on='gene_id', how='inner')

    # Extract contig ID
    df_all['Contig'] = df_all['gene_id'].apply(lambda x: '_'.join(x.split('_')[:-1]))

    # Harmonize column names
    df_all.rename(columns={'Sequence_coverage': 'Coverage'}, inplace=True)
    df_all['Identity'] = pd.NA
    df_all.drop(columns=['#Hit_Id', 'gene_id'], inplace=True)

    # Ensure Rep has same columns
    df_rep = df_rep[df_all.columns] if not df_rep.empty else pd.DataFrame(columns=df_all.columns)

    # Concatenate all annotations
    df_all = pd.concat([df_all, df_rep], ignore_index=True)

    # Final cleanup of 'tag' column
    df_all.rename(columns={'Reference_system': 'tag'}, inplace=True)

    # Write to output file
    df_all[['Contig', 'start', 'end', 'strand', 'Gene', 'Identity', 'Coverage', 'i-evalue', 'tag']].to_csv(output_path, sep="\t", index=False)
    print(f"Biomarker table written to: {output_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Retrieve plasmid biomarkers from COPLA output')
    parser.add_argument('-p', '--gene_pos', required=True, help='Path to gene_positions.tsv')
    parser.add_argument('-m', '--mob', required=True, help='Path to MOBScan output')
    parser.add_argument('-r', '--rep', required=True, help='Path to Pfinder output')
    parser.add_argument('-c', '--conj', default=None, help='Path to CONJScan output')
    parser.add_argument('-o', '--output', required=True, help='Output biomarkers path')
    args = parser.parse_args()

    main(args.gene_pos, args.mob, args.rep, args.conj, args.output)
