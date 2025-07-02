#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import argparse
import logging
from pathlib import Path
from typing import Tuple, List

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='prophage_parser.py is part of PITISfinder.')
    parser.add_argument('-p', '--provirus', required=True, type=Path)
    parser.add_argument('-t', '--taxonomy', required=True, type=Path)
    parser.add_argument('-a', '--ann_file', required=True, type=Path)
    parser.add_argument('-o', '--out_dir', required=True, type=Path)
    return parser.parse_args()


def extract_region(row: pd.Series, input_gbk: Path, out_dir: Path) -> Tuple[List[str], str, List[str]]:
    outname = row['Name']
    contig_id = row['Contig']
    amr_list = []
    vf_list = []
    df_list = []
    new_features = []

    for record in SeqIO.parse(input_gbk, "genbank"):
        if record.id != contig_id:
            continue

        start = int(row['Start']) - 1
        end = int(row['End'])
        sub_seq = record.seq[start:end]

        for feature in record.features:
            if feature.location.end < start or feature.location.start > end:
                continue

            # Adjust feature location to local coordinates
            new_start = max(feature.location.start, start) - start
            new_end = min(feature.location.end, end) - start

            new_location = FeatureLocation(new_start, new_end, strand=feature.location.strand)
            new_qualifiers = dict(feature.qualifiers)

            if feature.type == "CDS" and "AMR" in feature.qualifiers.get("tag", [""]):
                amr_list.append(feature.qualifiers.get("gene", [""])[0])
            elif feature.type == "CDS" and "VF" in feature.qualifiers.get("tag", [""]):
                vf_list.append(feature.qualifiers.get("gene", [""])[0])
            elif feature.type == "CDS" and "DF" in feature.qualifiers.get("tag", [""]):
                df_list.append(feature.qualifiers.get("gene", [""])[0])

            new_features.append(SeqFeature(location=new_location, type=feature.type, qualifiers=new_qualifiers))

        new_record = SeqRecord(
            seq=sub_seq,
            id=outname,
            name=outname,
            description=f"Prophage {outname}. At position {start}:{end} in {record.id} from {row['Sample']}",
            annotations={"molecule_type": "DNA"},
            features=new_features
        )

        # Write GenBank and FASTA
        with open(out_dir / f"{outname}.gbk", "w") as gbk_out:
            SeqIO.write(new_record, gbk_out, "genbank")

        with open(out_dir / f"{outname}.fasta", "w") as fasta_out:
            SeqIO.write(new_record, fasta_out, "fasta")

        row = row.copy()
        row['AMR'] = ';'.join(amr_list)
        row['VF'] = ';'.join(vf_list)
        row['DF'] = ';'.join(df_list)
        return row, outname

    raise ValueError(f"Contig ID '{contig_id}' not found in {input_gbk}")


def reformat_tables(df_provirus, df_taxonomy):
    df_info = df_provirus.merge(df_taxonomy, on="seq_name")
    df_info.columns = [c.capitalize() for c in df_info.columns]
    df_info.rename(columns={'Source_seq': 'Contig'
                    }, inplace=True)

    df_info['Sample'] = sample
    df_info['Name'] = "phage_" + df_info['Taxid'].astype(str) + "_" + df_info.index.astype(str) + f"_{sample}"
    df_info['LastLineage'] = (
            df_info['Lineage']
            .str.split(';')
            .apply(lambda x: next((i for i in reversed(x) if i.strip()), ''))
        )
    return df_info[['Sample', 'Contig', 'Name', 'Length', 'Start', 'End', 'LastLineage']]

if __name__ == "__main__":
    args = get_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    sample = '_'.join(args.ann_file.stem.split('_')[:-2])

    args.out_dir.mkdir(parents=True, exist_ok=True)
    report_out = args.out_dir / "prophage_summary.tsv"

    try:
        df_provirus = pd.read_table(args.provirus, header=0)
        df_taxonomy = pd.read_table(args.taxonomy, header=0)
    except Exception as e:
        logging.warning(f"{sample}: cannot open input file. {e}")
        exit(1)

    df_info = reformat_tables(df_provirus, df_taxonomy)

    columns = ['Sample', 'Contig', 'Name', 'LastLineage', 'Length', 'Start', 'End', 'AMR', 'VF', 'DF']
    summary_records = []

    for _, row in df_info.iterrows():
        try:
            updated_row, outname = extract_region(row, args.ann_file, args.out_dir)
            summary_records.append(updated_row[columns])
            single_tsv_path = args.out_dir / f"{outname}.tsv"
            updated_row[columns].to_frame().T.to_csv(single_tsv_path, sep='\t', index=False)
        except Exception as e:
            logging.warning(f"{row['Sample']}:{row['Contig']} - Extraction failed: {e}")

    summary_df = pd.DataFrame(summary_records, columns=columns)
    summary_df.to_csv(report_out, sep='\t', index=False)
    logging.info(f"Report saved to: {report_out}")
