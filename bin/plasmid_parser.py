#!/usr/bin/env python3

"""
Script to parse plasmid information obtained from MOBrecon and COPLA, merged with user AMR annotation.
"""

import pandas as pd
from Bio import SeqIO
from BCBio import GFF
import re
import argparse
import logging
from pathlib import Path
from typing import Tuple, Dict, Any, Union, List

def get_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        prog = 'plasmid_parser.py',
        description = 'plasmid_parser.py is part of PITISfinder.'
    )

    parser.add_argument(
        '-m', '--mobsuite', required=True,
        type=Path,
        help="Required. mobtyper_results.txt file"
    )

    parser.add_argument(
        '-q', '--qry_info', required=True,
        type=Path,
        help="Required. Copla qry_info.tsv output file"
    )

    parser.add_argument(
        '-p', '--ptu', required=True,
        type=Path,
        help="Required. Copla ptu_prediction.tsv output file"
    )

    parser.add_argument(
        '-a', '--ann_file', required=True,
        type=Path,
        help="Required. Annotation (.gff, .gff3) file to be processed"
    )

    parser.add_argument(
        '-r', '--res_file', required=True,
        type=Path,
        help="Required. Resistance annotation file (.tab) to be processed"
    )

    parser.add_argument(
        '-n', '--plasmid_name', required=True,
        type=str,
        help="Required. Plasmid name"
    )

    parser.add_argument(
        '-o', '--out_dir', required=True,
        type=Path,
        help='Required. Final output folder for reports'
    )

    return parser.parse_args()

def load_tsv(path: Path, **kwargs) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep='\t', **kwargs)
    except Exception as e:
        logging.error(f"Failed to read {path}: {e}")
        raise

def extract_plasmid_row(df: pd.DataFrame, col: str, value: str) -> pd.Series:
    parts = df[col].str.split(':', n=1, expand=True)
    df[['sample', 'id']] = parts
    filtered = df[df['id'] == value]
    if filtered.empty:
        raise ValueError(f"No records found for plasmid id '{value}' in {col}")
    return filtered.iloc[0]

def build_info(mrow: pd.Series, qrow: pd.Series, prow: pd.Series) -> Dict[str, any]:
    return {
        'mobsuite_id': mrow['id'],
        'ptu': prow.get('#Predicted_PTU'),
        'ptu_score': prow.get('Score'),
        'size': mrow.get('size'),
        'rep_mobsuite': mrow.get('rep_type(s)'),
        'rep_copla': qrow.get('Replicon'),
        'relaxase_mobsuite': mrow.get('relaxase_type(s)'),
        'mob_copla': qrow.get('MOB'),
        'mpf_mobsuite': mrow.get('mpf_type'),
        'mpf_copla': qrow.get('MPF'),
        'orit_mobsuite': mrow.get('orit_type(s)'),
        'predicted_mobility': mrow.get('predicted_mobility'),
        'primary_cluster_id': mrow.get('primary_cluster_id'),
        'secondary_cluster_id': mrow.get('secondary_cluster_id'),
        'mash_neighbor_identification': mrow.get('mash_neighbor_identification'),
        'mash_neighbor_distance': mrow.get('mash_neighbor_distance'),
        'AMR': qrow.get('AMR')
    }

def main():
    args = get_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


    # Output
    args.out_dir.mkdir(parents=True, exist_ok=True)
    report_out = args.out_dir / "plasmids_summary.tsv"
    pl_id = args.plasmid_name.split('_')[0]

    # Load tables
    df_mobs = load_tsv(args.mobsuite)
    df_qry = load_tsv(args.qry_info)
    df_ptu = load_tsv(args.ptu)
    # df_res = load_tsv(args.res_file)

    cols_keep = ['sample_id', 'size', 'rep_type(s)', 'relaxase_type(s)',
                 'mpf_type', 'orit_type(s)', 'predicted_mobility',
                 'primary_cluster_id', 'secondary_cluster_id',
                 'mash_neighbor_distance', 'mash_neighbor_identification']
    df_mobs = df_mobs[cols_keep]

    mrow = extract_plasmid_row(df_mobs, 'sample_id', pl_id)
    qrow = df_qry.iloc[0]
    prow = df_ptu.iloc[0]

    info = build_info(mrow, qrow, prow)
    df_report = pd.DataFrame([info])

    # Save report
    df_report.to_csv(report_out, sep='\t', index=False)
    logging.info(f"Report saved to: {report_out}")

if __name__ == "__main__":
    main()