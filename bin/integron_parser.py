#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re
import argparse
import logging
from pathlib import Path
from typing import Tuple, Dict, Any, List
from collections import defaultdict

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='integron_parser.py is part of PITISfinder.')
    parser.add_argument('-i', '--integron_file', required=True, type=Path)
    parser.add_argument('-a', '--ann_file', required=True, type=Path)
    parser.add_argument('--max_cas', type=int, default=12)
    parser.add_argument('-o', '--out_dir', required=True, type=Path)
    return parser.parse_args()

def basic_info(subdf: pd.DataFrame) -> Tuple[Dict[str, Any], pd.DataFrame]:
    integrase = subdf[subdf['annotation'].str.contains('intI')].iloc[0]
    if integrase['strand'] > 0:
        subdf = subdf[::-1]
    start, end = subdf['pos_beg'].min(), subdf['pos_end'].max()
    return {
        'type': integrase['type'],
        'integrase_model': integrase['model'],
        'strand': integrase['strand'],
        'start': start,
        'end': end,
        'size': end - start
    }, subdf

def merge_info(summary_df: pd.DataFrame, sample: str, contig: int, integron_name: str, d_info: Dict, cassettes: List, max_cassettes: int, amr_list: List, vf_list: List, df_list: List) -> pd.DataFrame:
    cassettes = (cassettes[:max_cassettes] + [None] * (max_cassettes - len(cassettes)))
    row = [sample, contig, integron_name, d_info['type'], d_info['size'], d_info['start'],
           d_info['end'], d_info['integrase_model']] + cassettes  + [','.join(amr_list)] + [','.join(vf_list)] + [','.join(df_list)]
    summary_df.loc[len(summary_df)] = row
    return row, summary_df

def extract_region(input_gbk: Path, contig_id: str, start: int, end: int, sample: str, count: int, out_dir: Path) -> List:
    for record in SeqIO.parse(input_gbk, "genbank"):
        if record.id != contig_id:
            continue

        new_features = []
        cds_names_list = []
        amr_list = []
        vf_list = []
        df_list = []

        new_start = start
        new_end = end
        for feature in record.features:

            if feature.location.start > end or feature.location.end < start:
                continue

            new_start = min(feature.location.start, new_start)
            new_end = max(feature.location.end, new_end)

            fstart = feature.location.start - new_start
            fend = feature.location.end - new_start
            new_loc = FeatureLocation(fstart, fend, strand=feature.location.strand)

            new_feature = SeqFeature(location=new_loc, type=feature.type, qualifiers=feature.qualifiers)
            new_features.append(new_feature)

            if feature.type == "CDS" and "AMR" in feature.qualifiers.get("tag", [""]):
                amr_list.append(feature.qualifiers.get("gene", [""])[0])

            elif feature.type == "CDS" and "VF" in feature.qualifiers.get("tag", [""]):
                vf_list.append(feature.qualifiers.get("gene", [""])[0])

            elif feature.type == "CDS" and "DF" in feature.qualifiers.get("tag", [""]):
                df_list.append(feature.qualifiers.get("gene", [""])[0])

            elif feature.type == "CDS" and "inti" not in feature.qualifiers.get('tag', [""])[0].lower():
                cds_names_list.append(feature.qualifiers.get("gene", ["protein"])[0])

        cassettes = '_'.join([re.sub(r'[^a-zA-Z0-9]', '', i) for i in cds_names_list])

        outname = f"int_{cassettes}_{sample}_{count}"
        new_record = SeqRecord(
            seq=record.seq[new_start:new_end],
            id=outname,
            annotations={"molecule_type": "DNA"},
            name=record.name,
            description=f"Integron {outname}. At position {new_start}:{new_end} in {record.id} from {sample}.",
            features=new_features
        )

        # Write GenBank and FASTA
        with open(out_dir / f"{outname}.gbk", "w") as gbk_out:
            SeqIO.write(new_record, gbk_out, "genbank")

        with open(out_dir / f"{outname}.fasta", "w") as fasta_out:
            SeqIO.write(new_record, fasta_out, "fasta")

        return cds_names_list, outname, amr_list, vf_list, df_list

    raise ValueError(f"Contig ID '{contig_id}' not found in {input_gbk}")

if __name__ == "__main__":
    args = get_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    sample = args.integron_file.stem

    args.out_dir.mkdir(parents=True, exist_ok=True)
    report_out = args.out_dir / "integron_summary.tsv"

    try:
        df_integrons = pd.read_table(args.integron_file, comment='#')
    except Exception as e:
        logging.warning(f"{sample}: cannot open integron file. {e}")
        exit()

    df_integrons = df_integrons[df_integrons['type'] != 'CALIN']
    if df_integrons.empty:
        logging.info(f"{sample}: no complete integrons found.")
        exit()

    columns = ['Sample', 'Contig', 'Name', 'Type', 'Length', 'Start', 'End', 'Integrase'] + [f'Cassette {i+1}' for i in range(args.max_cas)] + ['AMR', 'VF', 'DF']
    summary_df = pd.DataFrame(columns=columns)

    df_grouped = df_integrons.groupby(['ID_replicon', 'ID_integron'])

    for count, ((contig, integron), subdf) in enumerate(df_grouped):
        d_info, subdf = basic_info(subdf)
        if d_info['type'] == 'complete':
            cassettes, integron_name, amr_list, vf_list, df_list = extract_region(args.ann_file, contig, d_info['start'], d_info['end'], sample, count, args.out_dir)
            row, summary_df = merge_info(summary_df, sample, contig, integron_name, d_info, cassettes, args.max_cas, amr_list, vf_list, df_list)
            int_df = pd.DataFrame([row], columns=columns)
            int_df.to_csv(f'{args.out_dir / integron_name}.tsv', index=False, sep='\t')
        else:
            cassettes = []
            integron_name = f'in0_{sample}_{count}'
            row, summary_df = merge_info(summary_df, sample, contig, integron_name, d_info, cassettes, args.max_cas, [], [], [])

    summary_df.to_csv(report_out, sep='\t', index=False)
    logging.info(f"Report saved to: {report_out}")
