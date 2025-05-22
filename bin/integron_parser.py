#!/usr/bin/env python3

"""
Script to parse integrons obtained with Integron_parser, including protein and AMR annotations.
"""

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re
import argparse
import logging
from pathlib import Path
from typing import Tuple, Dict, Any, List

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
        '-a', '--ann_file', required=True,
        type=Path,
        help="Required. Annotation (.gff, .gff3) file to be processed"
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
        'type': integrase['type'],
        'integrase_model': integrase['model'],
        'strand': integrase['strand'],
        'start': start,
        'end': end,
        'size': end - start
    }, subdf

def merge_info(summary_df: pd.DataFrame, sample: str, contig: int, integron_name: str, d_info: Dict, cassettes: List, max_cassettes: int, amr_list: List) -> pd.DataFrame:
    cassettes = (cassettes[:max_cassettes] + [None] *
                (max_cassettes - len(cassettes)))                             # Truncate cassettes if too many, add None if too few
    row = [sample, contig, integron_name, d_info['type'], d_info['size'], d_info['start'],
           d_info['end'], d_info['integrase_model']] + cassettes  + [','.join(amr_list)]  # Create row to append
    summary_df.loc[len(summary_df)] = row
    return row, summary_df


def extract_region_with_attC(input_gbk: Path, contig_id: str, start: int, end: int, sample: str, count: int, attc_df: pd.DataFrame, out_dir: Path) -> List:
    for record in SeqIO.parse(input_gbk, "genbank"):
        if record.id != contig_id:
            continue

        sub_seq = record.seq[start:end]

        # Adjust and retain features within the region
        new_features = []
        cds_names_list = []
        amr_list = []
        for feature in record.features:
            if feature.location.end > end + 50 or feature.location.start < start - 50:
                continue

            new_start = max(feature.location.start, start) - start
            new_end = min(feature.location.end, end) - start
            new_loc = FeatureLocation(new_start, new_end, strand=feature.location.strand)

            new_feature = SeqFeature(location=new_loc, type=feature.type, qualifiers=feature.qualifiers)
            new_features.append(new_feature)

            if feature.type == "CDS" and "AMR" in feature.qualifiers.get("tag", [""]):
                gene_name = feature.qualifiers.get("gene", [""])[0]
                amr_list.append(gene_name)

            if feature.type == "CDS" and "inti" not in feature.qualifiers.get("gene", [""])[0].lower():
                gene_name = feature.qualifiers.get("gene", ["protein"])[0]
                cds_names_list.append(gene_name)
        cassettes = '_'.join([re.sub(r'[^a-zA-Z0-9]', '', i) for i in cds_names_list])

        # Add attC features from external table
        attc_subset = attc_df[
            (attc_df["ID_replicon"] == contig_id) &
            (attc_df["type_elt"].str.lower() == "attc") &
            (attc_df["pos_end"] > start) &
            (attc_df["pos_beg"] < end)
        ]
        for _, row in attc_subset.iterrows():
            attc_start = int(row["pos_beg"])
            attc_end = int(row["pos_end"])
            strand = row["strand"]
            model = row.get("model", "NA")

            local_start = max(attc_start, start) - start
            local_end = min(attc_end, end) - start

            location = FeatureLocation(local_start, local_end, strand=strand)
            attc_feature = SeqFeature(location=location, type="regulatory", qualifiers={"note": ["attC"], "model": [model]})
            new_features.append(attc_feature)

        # Create new record
        outname = f"int_{cassettes}_{sample}_{count}"
        new_record = SeqRecord(
            seq=sub_seq,
            id=outname,
            annotations={"molecule_type": "DNA"},
            name=record.name,
            description=f"Integron {outname}. At position {start}:{end} in {record.id} from {sample}.",
            features=new_features
        )

        # Sort features by ascending order
        source_feats = [f for f in new_record.features if f.type == 'source']
        other_feats = [f for f in new_record.features if f.type != 'source']
        other_feats.sort(key=lambda f: int(f.location.start))
        new_record.features = source_feats + other_feats

        # Write GenBank
        with open(out_dir / f"{outname}.gbk", "w") as gbk_out:
            SeqIO.write(new_record, gbk_out, "genbank")

        # Write FASTA
        with open(out_dir / f"{outname}.fasta", "w") as fasta_out:
            SeqIO.write(new_record, fasta_out, "fasta")

        return cds_names_list, outname, amr_list

    raise ValueError(f"Contig ID '{contig_id}' not found in {input_gbk}")

if __name__ == "__main__":
    args = get_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    sample = args.integron_file.stem

    # Output
    args.out_dir.mkdir(parents=True, exist_ok=True)
    report_out = args.out_dir / "integron_summary.tsv"

    # Try open integron file, quit if fail
    try:
        df_integrons = pd.read_table(args.integron_file, comment='#')
    except Exception as e:
        logging.warning(f"{sample}: cannot open integron file. {e}")
        exit()

    # Keep only complete integrons, quit if none
    df_integrons = df_integrons[df_integrons.type != 'CALIN']
    if df_integrons.empty:
        logging.info(f"{sample}: no complete integrons found.")
        exit()

    # Empty dataframe with {max_cassettes} cassettes
    columns = ['Sample', 'Contig', 'Name', 'Type', 'Size', 'Start', 'End', 'Integrase'] + [f'Cassette {i+1}' for i in range(args.max_cas)] + ['AMR']
    summary_df = pd.DataFrame(columns=columns)

    # Group by integron number and contig
    df_grouped = df_integrons.groupby(['ID_replicon', 'ID_integron'])

    for count, ((contig, integron), subdf) in enumerate(df_grouped):    # keep count of integron number
        d_info, subdf = basic_info(subdf)                               # Read coordinates
        if d_info['type'] == 'complete':
            cassettes, integron_name, amr_list = extract_region_with_attC(args.ann_file, contig, d_info['start'], d_info['end'], sample, count, subdf, args.out_dir)
            row, summary_df = merge_info(summary_df, sample, contig, integron_name, d_info, cassettes, args.max_cas, amr_list)
            int_df = pd.DataFrame([row], columns=columns)
            int_df.to_csv(f'{args.out_dir / integron_name}.tsv', index=False, sep='\t')
        else:
            cassettes = []
            integron_name = f'in0_{sample}_{count}'
            row, summary_df = merge_info(summary_df, sample, contig, integron_name, d_info, cassettes, args.max_cas, [])

    summary_df.to_csv(report_out, sep='\t', index=False)
    logging.info(f"Report saved to: {report_out}")
