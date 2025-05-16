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
        '-m', '--mobsuite_typer', required=True,
        type=Path,
        help="Required. mobtyper_results.txt file"
    )

    parser.add_argument(
        '-r', '--mobsuite_report', required=True,
        type=Path,
        help="Required. contig_report.txt file"
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
        '-g', '--gbk', required=True,
        type=Path,
        help="Required. Annotation (.gbk, .gbff) file to be processed"
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

def extract_plasmid_row(df: pd.DataFrame, value: str) -> pd.Series:
    filtered = df[df['primary_cluster_id'] == value]
    return filtered

def build_info(mrow: pd.Series, contigs:List, qrow: pd.Series, prow: pd.Series, amr_genes: List) -> Dict[str, any]:
    return {
        'sample': mrow['sample_id'],
        'contig': ','.join(contigs),
        'mobsuite_id': mrow['primary_cluster_id'],
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
        'AMR': ','.join(amr_genes)
    }

def main():
    args = get_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Output
    args.out_dir.mkdir(parents=True, exist_ok=True)
    report_out = args.out_dir / f"{args.plasmid_name}.tsv"
    gbk_out = args.out_dir / f"{args.plasmid_name}.gbk"
    pl_id = args.plasmid_name.split('_')[0]

    # Load tables
    df_mobt = load_tsv(args.mobsuite_typer)
    parts = df_mobt['sample_id'].str.split(':', n=1, expand=True)
    df_mobt[['sample_id', 'primary_cluster_id']] = parts
    df_mobr = load_tsv(args.mobsuite_report)
    df_qry = load_tsv(args.qry_info)
    df_ptu = load_tsv(args.ptu)
    # df_res = load_tsv(args.res_file)

    mobt_filt = extract_plasmid_row(df_mobt, pl_id)
    print(df_mobr)
    mobr_filt = extract_plasmid_row(df_mobr, pl_id)
    print(mobr_filt)
    contigs = mobr_filt['contig_id'].tolist()

    qrow = df_qry.iloc[0]
    prow = df_ptu.iloc[0]

    # Save annotation (GBK)
    selected_records = []
    with args.gbk.open("r") as gbk_file:
        for record in SeqIO.parse(gbk_file, "genbank"):
            if record.id in contigs:
                selected_records.append(record)

    with gbk_out.open("w") as out_handle:
        SeqIO.write(selected_records, out_handle, "genbank")
    logging.info(f"Filtered GenBank file saved to: {gbk_out}")

    # Retrieve AMR genes from GBK annotation
    amr_genes = []
    for record in selected_records:
        for feature in record.features:
            if feature.type == "CDS" and 'AMR' in feature.qualifiers.get('tag', [""]):
                gene_name = feature.qualifiers.get("gene", [""])[0]
                amr_genes.append(gene_name)

    # Retrieve report info
    info = build_info(mobt_filt.iloc[0], contigs, qrow, prow, amr_genes)
    df_report = pd.DataFrame([info])

    # Save report
    df_report.to_csv(report_out, sep='\t', index=False)
    logging.info(f"Report saved to: {report_out}")


if __name__ == "__main__":
    main()