#!/usr/bin/env python3

"""
Script to merge GBK and TAB (from Integron_finder) biomarkers.
Adds and modifies features to facilitate parsing and annotation.
"""

import argparse
import pandas as pd
import numpy as np
import logging
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def load_tab(tab_path):
    df = pd.read_csv(tab_path, sep="\t", header=0, comment="#")
    df['pos_beg'] = df['pos_beg'].astype(int)
    df['pos_end'] = df['pos_end'].astype(int)
    df['strand'] = df['strand'].astype(int)
    return df

def reformat_table(df_biom):
    index_protein = df_biom[df_biom['annotation'] == "protein"].index
    df_biom.drop(index_protein, inplace=True)
    df_biom.sort_values(by=['ID_replicon', 'pos_beg'], inplace=True)
    df_biom['tag'] = df_biom['annotation'].str.split('_').str[0]
    df_biom['product'] = df_biom['tag']
    df_biom['gene'] = df_biom['tag']
    d_gbtype = {
        "attC": "regulatory",
        "attI": "regulatory",
        "Promoter": "regulatory",
        "protein": "CDS"
    }
    df_biom['gb_type'] = df_biom['type_elt'].map(d_gbtype)
    return df_biom

def annotate_record(record, df, nts_diff):
    # select only tab hits for this contig
    df_rec = df[df['ID_replicon'] == record.id]
    if df_rec.empty:
        return record

    # Remove any existing features that overlap any tab region
    new_features = []
    for feat in record.features:
        fstart = int(feat.location.start) + 1 + nts_diff
        fend   = int(feat.location.end) - nts_diff
        # check overlap with any tab row
        overlap = False
        for i, row in df_rec.iterrows():
            if not (fend < row['pos_beg'] or fstart > row['pos_end']):
                overlap = True
                df_rec.loc[i, 'product'] = feat.qualifiers.get('product', [row['product']])[0]
                df_rec.loc[i, 'gene'] = feat.qualifiers.get('gene', [row['gene']])[0]
                break
        if not overlap:
            new_features.append(feat)
    record.features = new_features

    # Add/overwrite features from tab (prioritized)
    for _, row in df_rec.iterrows():
        loc = FeatureLocation(row['pos_beg'] - 1, row['pos_end'], strand=row['strand'])
        qualifiers = {
            'db_xref': [f"Integron_finder:{row['annotation']}_{row['model']}"],
            'product': [f"{row['product']}, {row['model']} model"],
            'locus_tag': [f"{row['ID_replicon']}_{row['pos_beg']}_{row['pos_end']}"],
            'protein_id': [f"gnl|Integron_finder|{row['ID_replicon']}_{row['pos_beg']}_{row['pos_end']}"],
            'tag': row['tag'],
            'codon_start': ['1'],
            'transl_table': ['11'],
            'inference': [f"Integron_finder prediction, {row['product']}, {row['model']}"],
            'gene': [row['gene']],
            'mge_element': ['yes']
        }
        new_feat = SeqFeature(location=loc, type=row['gb_type'], qualifiers=qualifiers)
        record.features.append(new_feat)

    # reorder features: keep 'source' first, then CDS (and other) sorted by start
    source_feats = [f for f in record.features if f.type == 'source']
    other_feats = [f for f in record.features if f.type != 'source']
    other_feats.sort(key=lambda f: int(f.location.start))
    record.features = source_feats + other_feats
    return record

def main(biomark, gbk, output, nts_diff):
    try:
        df_biom = load_tab(biomark)
    except Exception as e:
        logging.warning(f"No integrons found in sample. {e}")
        exit(0)

    df_annotation = reformat_table(df_biom)
    records = list(SeqIO.parse(gbk, 'genbank'))
    annotated = []
    for rec in records:
        annotated.append(annotate_record(rec, df_annotation, nts_diff))
    SeqIO.write(annotated, output, 'genbank')
    print(f"Wrote annotated GenBank to {output}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge Integron_finder biomarkers with GenBank file annotations')
    parser.add_argument('-g', '--gbk', required=True, help='Path to input GenBank file')
    parser.add_argument('-b', '--biomarkers', required=True, help='Path to sample.integrons file')
    parser.add_argument('-o', '--output', default='annotated_output.gbk', help='Output GenBank path')
    parser.add_argument('-n', '--nts_diff', default=15, help='Allowed N nucleotide overlapping (default: 15)')
    args = parser.parse_args()

    main(args.biomarkers, args.gbk, args.output, args.nts_diff)
