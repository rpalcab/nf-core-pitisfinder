#!/usr/bin/env python3

"""
Script to solve overlapping BLAST hits from a .tsv file.
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import logging

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog = 'blastn_overlap_solver.py',
        description = 'Solves overlapping hits from Blastn search.'
    )

    parser.add_argument(
        '-i', '--input_file', required=True,
        type=Path,
        help="Required. Input file to be filtered"
    )
    parser.add_argument(
        '-o', '--out', required=True,
        type=Path,
        help='Required. Output directory.'
    )

    return parser.parse_args()

def process_table(df: pd.DataFrame) -> pd.DataFrame:

    header = [
        'IS', 'contig', 'qstart', 'qend', 'qlen',
        'sstart0', 'send0', 'slen', 'pident', 'qcovhsp',
        'length', 'mismatch', 'score', 'evalue'
        ]
    df.columns = header
    # Contigs as numbers
    # df['contig'] = df['contig'].str.replace(r'^[a-zA-Z]+_', '', regex=True).astype(int)
    # Redefines sequence coordinates in negative frames
    coords = np.sort(df[['sstart0', 'send0']], axis=1)
    df['sstart'], df['send'] = coords[:, 0], coords[:, 1]
    df.drop(columns=['sstart0', 'send0'], inplace=True)

    return df

def remove_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    grouped_filtered = []

    for query, subdf in df.groupby("IS", sort=False):
        subdf = subdf.sort_values(
            by=["send", "score"],
            ascending=[True, False]
        )
        subdf = subdf.drop_duplicates(
            subset=["send"],
            keep="first"
        )
        subdf = subdf.sort_values(
            by=["sstart", "score"],
            ascending=[True, False]
        )
        subdf = subdf.drop_duplicates(
            subset=["sstart"],
            keep="first"
        )
        grouped_filtered.append(subdf)

    return pd.concat(grouped_filtered).reset_index(drop=True)

def resolve_overlaps(df: pd.DataFrame) -> pd.DataFrame:
    filtered = []
    for _, row in df.iterrows():
        # Checks hits are in the same IS and the current hit starts before the previous one ends
        while filtered and row.IS == filtered[-1].IS and row.sstart < filtered[-1].send:
            if row.score > filtered[-1].score:
                filtered.pop()
            else:
                break
        else:
            filtered.append(row)

    return pd.DataFrame(filtered)

def main():
    args = get_args()

    # Prepare
    args.out.mkdir(parents=True, exist_ok=True)
    out_file = args.out / f"{args.input_file.stem}_filtered.tsv"
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    # Input processing
    df = pd.read_table(args.input_file, sep='\t')
    logging.info(f"Initial number of hits: {len(df)}")
    df = process_table(df)
    df = remove_duplicates(df)
    filtered_df = resolve_overlaps(df)

    # Output
    filtered_df.to_csv(out_file, sep='\t', index=False)
    logging.info(f"Wrote {len(filtered_df)} resolved hits to: {out_file}")

if __name__ == '__main__':
    main()
