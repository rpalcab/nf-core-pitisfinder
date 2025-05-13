#!/usr/bin/env python3

"""
Script to parse integrons obtained with Integron_parser, including protein and AMR annotations.
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
        prog = 'integron_parser.py',
        description = 'integron_parser.py is part of PITISfinder.'
    )

    parser.add_argument(
        '-i', '--integron_file', required=True,
        type=Path,
        help="Required. Integron file to be filtered and processed"
    )

    parser.add_argument(
        '-f', '--fasta_file', required=True,
        type=Path,
        help="Required. Fasta file to be processed"
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
        '-s', '--sample', required=True,
        type=str,
        help="Required. Sample name"
    )

    parser.add_argument(
        '--nts_diff', type=int, default=10,
        help="Coordinates differences tolerated when merging annotations. Default: 10"
    )

    parser.add_argument(
        '--id_perc', type=int, default=90,
        help="Minimum ID percentage required for antimicrobial resistance identification. Default: 90"
    )

    parser.add_argument(
        '--q_cov', type=int, default=80,
        help="Minimum query coverage required for antimicrobial resistance identification. Default: 80"
    )

    parser.add_argument(
        '--max_cas', type=int, default=12,
        help="Maximum cassettes in integron. Default: 12"
    )

    parser.add_argument(
        '-o', '--out_dir', required=True,
        type=Path,
        help='Required. Final output folder for reports'
    )

    return parser.parse_args()

def basic_info(subdf: pd.DataFrame) -> Tuple[Dict[str, Any], pd.DataFrame]:
    integrase = subdf[subdf['annotation'] == 'intI'].iloc[0]
    if integrase['strand'] > 0:
        subdf = subdf[::-1]

    start, end = subdf['pos_beg'].min(), subdf['pos_end'].max()
    return {
        'integrase_model': integrase['model'],
        'strand': integrase['strand'],
        'start': start,
        'end': end,
        'size': end - start
    }, subdf

def read_gff(ann_file: Path, contig: str, d_info: Dict) -> pd.DataFrame:
    with ann_file.open() as infile:
        records = list(GFF.parse(infile))

    data = [
        [rec.id, feat.type, feat.location.start, feat.location.end,
         feat.location.strand, feat.qualifiers.get('Name', [None])[0],
         feat.qualifiers.get('gene', [None])[0]]
        for rec in records for feat in rec.features
    ]

    df = pd.DataFrame(data, columns=[
        "seqid", "type", "start", "end",
        "strand", "name", "gene"
    ])
    return df[
        (df['type'] == 'CDS') &
        (df['seqid'] == contig) &
        (df['start'] > d_info['start']) &
        (df['end'] < d_info['end'])
    ]

def read_abr(res_file: Path, contig: str, d_info: Dict, pid: int, qcov:int) -> pd.DataFrame:
    df = pd.read_table(res_file, dtype={'START': int, 'END': int})
    df.drop(df.loc[
        (df['%COVERAGE'] < qcov) |
        (df['%IDENTITY'] < pid)
    ].index, inplace=True)     # Remove low quality hits
    df = df[
        (df.SEQUENCE == contig) &
        (df.START > d_info['start']) &
        (df.END < d_info['end'])
    ]  # Keep genes inside integron
    return df

def read_annotation(row: pd.Series, df: pd.DataFrame, start: str, end: str, annotation: str, nts_diff: int) -> Any:
    df_filtered = df[(abs(df[start] - row['pos_beg']) <= nts_diff) &    # Takes annotation in region
                     (abs(df[end] - row['pos_end']) <= nts_diff)]
    return df_filtered[annotation].iloc[0] if not df_filtered.empty else None   # Returns annotation if it exists

def update_annotation(row: pd.Series, df_ann: pd.DataFrame, df_abr: pd.DataFrame, nts_diff: int) -> Union[pd.DataFrame, pd.Series]:
    if row.annotation in {'attC', 'intI'}:
        return row.annotation

    return (
        read_annotation(row, df_abr, 'START', 'END', 'GENE', nts_diff) or
        read_annotation(row, df_ann, 'start', 'end', 'gene', nts_diff) or
        row.annotation
    )

def merge_annotations(subdf: pd.DataFrame, df_ann: pd.DataFrame, df_abr: pd.DataFrame, nts_diff: int) -> pd.DataFrame:
    subdf = subdf.copy()
    subdf['merged_annotation'] = subdf.apply(
        lambda row:update_annotation(row, df_ann, df_abr, nts_diff), axis=1
    )
    return subdf

def get_cassettes(unified_df: pd.DataFrame, count:int, sample: str) -> Tuple[List, str]:
    l_cassettes = unified_df.loc[(unified_df['type_elt'] == 'protein') &
                                 (unified_df['annotation'] != 'intI')]['merged_annotation'].to_list()
    cassettes = '_'.join([re.sub(r'[^a-zA-Z0-9]', '', i) for i in l_cassettes])
    name = f'int_{cassettes}_{sample}_{count}'
    return l_cassettes, name

def merge_info(summary_df: pd.DataFrame, sample: str, contig: int, integron_name: str, d_info: Dict, cassettes: List, max_cassettes: int) -> pd.DataFrame:
    cassettes = (cassettes[:max_cassettes] + [None] *
                (max_cassettes - len(cassettes)))                             # Truncate cassettes if too many, add None if too few
    row = [sample, contig, integron_name, d_info['size'], d_info['start'],
           d_info['end'], d_info['integrase_model']] + cassettes    # Create row to append
    summary_df.loc[len(summary_df)] = row
    return summary_df

def write_fasta(cds_output_file: Path, d_int: Dict) -> None:
    with open(cds_output_file, "w") as cds_output_handle:
        for k, v in d_int.items():
            cds_output_handle.write(f'>{k}' + "\n")
            cds_output_handle.write(v + "\n")
    return None

def save_sequence(unified_df: pd.DataFrame, fa_file: Path, out_dir: Path) -> None:
    cont = unified_df['Contig']
    start = unified_df['Start'] - 1
    end = unified_df['End'] - 1
    name = unified_df['Name']

    for seq_record in SeqIO.parse(fa_file, "fasta"):
        if seq_record.id == cont:
            d = {name: str(seq_record.seq[start:end])}
            write_fasta(out_dir / f'{name}.fasta', d)
            break
    return None



if __name__ == "__main__":
    args = get_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Output
    args.out_dir.mkdir(parents=True, exist_ok=True)
    report_out = args.out_dir / "integrons_summary.tsv"

    # Try open integron file, quit if fail
    try:
        df_integrons = pd.read_table(args.integron_file, comment='#')
    except Exception as e:
        logging.warning(f"{args.sample}: cannot open integron file. {e}")
        exit()

    # Keep only complete integrons, quit if none
    df_integrons = df_integrons[df_integrons.type == 'complete']
    if df_integrons.empty:
        logging.info(f"{args.sample}: no complete integrons found.")
        exit()

    # Empty dataframe with {max_cassettes} cassettes
    summary_df = pd.DataFrame(columns=['Sample', 'Contig', 'Name', 'Size', 'Start', 'End', 'Integrase'] +
                                      [f'Cassette {i+1}' for i in range(args.max_cas)])

    # Divide into integrons and chromosomes
    df_grouped = df_integrons.groupby(['ID_replicon', 'ID_integron'])

    for count, ((contig, integron), subdf) in enumerate(df_grouped):    # keep count of integron number
        d_info, subdf = basic_info(subdf)                               # Read coordinates
        df_ann = read_gff(args.ann_file, contig, d_info)                     # Read bakta annotation
        df_abr = read_abr(args.res_file, contig, d_info, args.id_perc, args.q_cov)   # Read abricate annotation
        unified_df = merge_annotations(subdf, df_ann, df_abr, args.nts_diff)    # Update annotations: abr > bakta > i_f
        cassettes, integron_name = get_cassettes(unified_df, count, args.sample)                     # get cassette, get name
        summary_df = merge_info(summary_df, args.sample, contig, integron_name, d_info, cassettes, args.max_cas)
        save_sequence(summary_df.iloc[-1], args.fasta_file, args.out_dir)                                    # Save integron sequence (read fasta, get coordinates, save sequence)

    summary_df.to_csv(report_out, index=False)
    logging.info(f"Report saved to: {report_out}")
