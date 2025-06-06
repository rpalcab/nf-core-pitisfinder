#!/usr/bin/env python3

"""
Script to merge GBK and TAB (from MOBrecon) biomarkers.
Adds and modifies features to facilitate parsing and annotation.
"""

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def load_tab(tab_path):
    df = pd.read_csv(tab_path, sep="\t", header=0)
    df['sstart'] = df['sstart'].astype(int)
    df['send'] = df['send'].astype(int)
    return df

def reformat_table(df_biom):
    df_biom.sort_values(by=['sseqid', 'sstart'], inplace=True)
    d_tag = {
        "mate-pair-formation": "MPF",
        "oriT": "oriT",
        "relaxase": "MOB",
        "replicon": "Replicon"
    }
    df_biom['tag'] = df_biom['biomarker'].map(d_tag)
    df_biom['product'] = df_biom['qseqid']
    df_biom['gene'] = df_biom['qseqid'].str.split('|').str[0]
    return df_biom

def annotate_record(record, df, nts_diff):
    # select only tab hits for this contig (by SEQUENCE)
    df_rec = df[df['sseqid'] == record.id]
    if df_rec.empty:
        return record

    # Remove any existing CDS features that overlap any tab region
    new_features = []
    for feat in record.features:
        if feat.type != 'CDS':
            new_features.append(feat)
            continue
        # existing CDS coords
        fstart = int(feat.location.start) + 1 + nts_diff
        fend   = int(feat.location.end) - nts_diff
        # check overlap with any tab row
        overlap = False
        for i, row in df_rec.iterrows():
            if not (fend < row['sstart'] or fstart > row['send']):
                overlap = True
                df_rec.loc[i, 'product'] = feat.qualifiers.get('product', [row['product']])[0]
                df_rec.loc[i, 'gene'] = feat.qualifiers.get('gene', [row['gene']])[0]
                break
        if not overlap:
            new_features.append(feat)
    record.features = new_features

    # Add/overwrite features from tab (prioritized)
    for _, row in df_rec.iterrows():
        strand_val = 1 if row['sstrand'] == 'plus' else -1
        loc = FeatureLocation(row['sstart'] - 1, row['send'], strand=strand_val)
        seq_segment = record.seq[loc.start:loc.end]
        if strand_val == -1:
            seq_segment = seq_segment.reverse_complement()
        translation = str(seq_segment.translate(table=11, to_stop=False))
        qualifiers = {
            'db_xref': [f"MOBsuite:{row['qseqid']}"],
            'product': [row['product']],
            'locus_tag': [f"{row['qseqid']}_{row['sstart']}_{row['send']}"],
            'protein_id': [f"gnl|MOBsuite|{row['qseqid']}_{row['sstart']}_{row['send']}"],
            'tag': row['tag'],
            'translation': [translation],
            'codon_start': ['1'],
            'transl_table': ['11'],
            'inference': [f"MOBsuite prediction, id% {row['pident']}, qcov% {row['qcovhsp']} to {row['gene']}"],
            'gene': [row['gene']],
            'mge_element': ['yes']
        }
        new_feat = SeqFeature(location=loc, type='CDS', qualifiers=qualifiers)
        record.features.append(new_feat)

    # reorder features: keep 'source' first, then CDS (and other) sorted by start
    source_feats = [f for f in record.features if f.type == 'source']
    other_feats = [f for f in record.features if f.type != 'source']
    other_feats.sort(key=lambda f: int(f.location.start))
    record.features = source_feats + other_feats
    return record

def main(biomark, gbk, output, nts_diff):
    df_biom = load_tab(biomark)
    df_annotation = reformat_table(df_biom)
    records = list(SeqIO.parse(gbk, 'genbank'))
    annotated = []
    for rec in records:
        annotated.append(annotate_record(rec, df_annotation, nts_diff))
    SeqIO.write(annotated, output, 'genbank')
    print(f"Wrote annotated GenBank to {output}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge MOBrecon biomarkers with GenBank file annotations')
    parser.add_argument('-g', '--gbk', required=True, help='Path to input GenBank file')
    parser.add_argument('-b', '--biomarkers', required=True, help='Path to biomarkers.blast.txt file')
    parser.add_argument('-o', '--output', default='annotated_output.gbk', help='Output GenBank path')
    parser.add_argument('-n', '--nts_diff', default=15, help='Allowed N nucleotide overlapping (default: 15)')
    args = parser.parse_args()

    main(args.biomarkers, args.gbk, args.output, args.nts_diff)
