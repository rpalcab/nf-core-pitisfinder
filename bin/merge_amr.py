#!/usr/bin/env python3

"""
Script to merge GBK and TAB (from Abricate) annotations.
Adds and modifies features to facilitate parsing and annotation.
"""

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def load_tab(tab_path):
    df = pd.read_csv(tab_path, sep="\t", header=0)
    df['START'] = df['START'].astype(int)
    df['END'] = df['END'].astype(int)
    return df

def annotate_record(record, df, nts_diff):
    # select only tab hits for this contig (by SEQUENCE)
    df_rec = df[df['SEQUENCE'] == record.id]
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
        for _, row in df_rec.iterrows():
            if not (fend < row['START'] or fstart > row['END']):
                overlap = True
                break
        if not overlap:
            new_features.append(feat)
    record.features = new_features

    # Add/overwrite features from tab (prioritized)
    for _, row in df_rec.iterrows():
        strand_val = 1 if row['STRAND'] == '+' else -1
        loc = FeatureLocation(row['START'] - 1, row['END'], strand=strand_val)
        seq_segment = record.seq[loc.start:loc.end]
        translation = str(seq_segment.translate(table=11, to_stop=True))
        qualifiers = {
            'db_xref': [f"{row['DATABASE']}:{row['ACCESSION']}]"],
            'product': [row['PRODUCT']],
            'locus_tag': [f"{row['ACCESSION']}_{row['START']}_{row['END']}"],
            'protein_id': [f"gnl|Abricate|{row['ACCESSION']}_{row['START']}_{row['END']}"],
            'tag': ['AMR'],
            'translation': [translation],
            'codon_start': ['1'],
            'transl_table': ['11'],
            'inference': [f"Abricate prediction, id% {row['%IDENTITY']}, qcov% {row['%COVERAGE']} to {row['ACCESSION']}"],
            'gene': [row['GENE']],
            'resistance': [row['RESISTANCE']]
        }
        new_feat = SeqFeature(location=loc, type='CDS', qualifiers=qualifiers)
        record.features.append(new_feat)

    # reorder features: keep 'source' first, then CDS (and other) sorted by start
    source_feats = [f for f in record.features if f.type == 'source']
    other_feats = [f for f in record.features if f.type != 'source']
    other_feats.sort(key=lambda f: int(f.location.start))
    record.features = source_feats + other_feats
    return record

def main(tab, gbk, output, nts_diff):
    df = load_tab(tab)
    records = list(SeqIO.parse(gbk, 'genbank'))
    annotated = []
    for rec in records:
        annotated.append(annotate_record(rec, df, nts_diff))
    SeqIO.write(annotated, output, 'genbank')
    print(f"Wrote annotated GenBank to {output}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge abricate.tab with GenBank file annotations')
    parser.add_argument('-t', '--tab', required=True, help='Path to abricate.tab file')
    parser.add_argument('-g', '--gbk', required=True, help='Path to input GenBank file')
    parser.add_argument('-o', '--output', default='annotated_output.gbk', help='Output GenBank path')
    parser.add_argument('-n', '--nts_diff', default=15, help='Allowed N nucleotide overlapping (default: 15)')
    args = parser.parse_args()
    main(args.tab, args.gbk, args.output, args.nts_diff)
