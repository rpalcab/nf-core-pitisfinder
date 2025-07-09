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
    df_biom['source'] = 'MOBsuite'
    return df_biom

def load_reformat_copla(copla):
    df_copla = pd.read_table(copla, sep='\t', header=0)
    adapted_df = pd.DataFrame({
        'qseqid': df_copla['Gene'],
        'sseqid': df_copla['Contig'],
        'qlen': pd.NA,
        'slen': pd.NA,
        'qstart': pd.NA,
        'qend': pd.NA,
        'sstart': df_copla['start'].astype(int),
        'send': df_copla['end'].astype(int),
        'length': pd.NA,
        'mismatch': pd.NA,
        'pident': df_copla['Identity'],
        'qcovhsp': df_copla['Coverage'],
        'qcovs': pd.NA,
        'sstrand': df_copla['strand'].map(lambda x: 'plus' if str(x) == '1' else 'minus' if str(x) == '-1' else pd.NA),
        'evalue': df_copla['i-evalue'],
        'bitscore': pd.NA,
        'biomarker': df_copla['tag'],
        'tag': df_copla['tag'],
        'product': df_copla['Gene'],
        'gene': df_copla['Gene'],
        'source': 'COPLA'
    })
    return adapted_df

def drop_overlapping_mobsuite(df_annotation, nts_diff=15):
    """
    Drops ovelapping hits (+- 15 nts difference).
    Prioritizes COPLA results over MOBsuite.
    """
    keep_rows = []
    mob_mask = df_annotation['source'] == 'MOBsuite'
    copla_mask = df_annotation['source'] == 'COPLA'
    df_mob_rows = df_annotation[mob_mask].copy()
    df_copla_rows = df_annotation[copla_mask].copy()

    for i, mob_row in df_mob_rows.iterrows():
        overlap = False
        for _, copla_row in df_copla_rows.iterrows():
            if mob_row['sseqid'] == copla_row['sseqid']:
                # Check overlap with nts_diff padding
                if not (mob_row['send'] + nts_diff < copla_row['sstart'] or
                        mob_row['sstart'] - nts_diff > copla_row['send']):
                    overlap = True
                    break
        if not overlap:
            keep_rows.append(i)

    # Keep all COPLA + non-overlapping MOBsuite. Drop overlapping COPLA hits
    df_concat = pd.concat([
        df_annotation[copla_mask],
        df_annotation.loc[keep_rows]
    ], ignore_index=True)
    df_concat.drop_duplicates(subset=['sseqid', 'sstart', 'send'], inplace=True)

    return df_concat

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
        # Infer strand for missing COPLA features overlapping this CDS
        fstart = int(feat.location.start) + 1 + nts_diff
        fend = int(feat.location.end) - nts_diff
        for i, row in df_rec[df_rec['source'] == 'COPLA'].iterrows():
            if row['sstrand'] not in ['plus', 'minus']:
                if not (fend < row['sstart'] or fstart > row['send']):
                    inferred_strand = feat.location.strand
                    df_rec.at[i, 'sstrand'] = 'plus' if inferred_strand == 1 else 'minus'
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
            'db_xref': [f"{row['source']}:{row['qseqid']}"],
            'product': [row['product']],
            'locus_tag': [f"{row['qseqid']}_{row['sstart']}_{row['send']}"],
            'protein_id': [f"gnl|{row['source']}|{row['qseqid']}_{row['sstart']}_{row['send']}"],
            'tag': row['tag'],
            'translation': [translation],
            'codon_start': ['1'],
            'transl_table': ['11'],
            'inference': [f"{row['source']} prediction, id% {row['pident']}, qcov% {row['qcovhsp']} to {row['gene']}"],
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

def main(mobsuite, copla, gbk, output, nts_diff):
    df_biom = load_tab(mobsuite)
    df_mobsuite = reformat_table(df_biom)
    df_copla = load_reformat_copla(copla)
    df_annotation = pd.concat([df_mobsuite, df_copla], ignore_index=True)
    df_annotation = drop_overlapping_mobsuite(df_annotation, nts_diff=int(nts_diff))
    records = list(SeqIO.parse(gbk, 'genbank'))
    annotated = []
    for rec in records:
        annotated.append(annotate_record(rec, df_annotation, nts_diff))
    SeqIO.write(annotated, output, 'genbank')
    print(f"Wrote annotated GenBank to {output}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge MOBrecon biomarkers with GenBank file annotations')
    parser.add_argument('-g', '--gbk', required=True, help='Path to input GenBank file')
    parser.add_argument('-m', '--mobsuite', required=True, help='Path to biomarkers.blast.txt file')
    parser.add_argument('-c', '--copla', required=True, help='Path to copla_biomarkers file')
    parser.add_argument('-o', '--output', default='annotated_output.gbk', help='Output GenBank path')
    parser.add_argument('-n', '--nts_diff', default=15, help='Allowed N nucleotide overlapping (default: 15)')
    args = parser.parse_args()

    main(args.mobsuite, args.copla, args.gbk, args.output, args.nts_diff)
