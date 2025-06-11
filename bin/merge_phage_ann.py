#!/usr/bin/env python3

"""
Script to merge GBK and TAB (from geNomad) biomarkers.
Adds and modifies features to facilitate parsing and annotation.
"""

import argparse
import pandas as pd
import numpy as np
import logging
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def load_tab(tab_path):
    df = pd.read_csv(tab_path, sep="\t", header=None, skiprows=[0])
    df.columns = ['gene', 'start', 'end', 'length', 'strand', 'gc_content', 'genetic_code',
                  'rbs_motif', 'marker', 'evalue', 'bitscore', 'uscg', 'taxid', 'taxname',
                  'annotation_accessions', 'annotation_description', 'unk1', 'unk2', 'domain', 'product']
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['strand'] = df['strand'].astype(int)
    df['contig'] = df['gene'].str.split('|').str[0]
    df['tag'] = "Phage protein"
    return df

def reformat_table(df_biom, df_int):
    integrases = set(df_int['integrases'])                                          # Annotate integrases
    mask = df_biom['gene'].isin(integrases)
    df_biom.loc[mask, 'annotation_accessions'] = -1
    df_biom.loc[mask, 'tag'] = 'Integrase'

    index_nophage = df_biom[df_biom['annotation_accessions'] == 1].index          # Remove everything not phage related
    df_biom.drop(index_nophage, inplace=True)

    df_biom = df_biom[~(df_biom['product'].isna() & (df_biom['tag'] == 'Phage protein'))]                        # Remove everything without product description

    df_biom.sort_values(by=['contig', 'start'], inplace=True)
    print(df_biom.head())
    return df_biom

def annotate_record(record, df, nts_diff):
    # select only tab hits for this contig
    df_rec = df[df['contig'] == record.id]
    if df_rec.empty:
        return record

    # Remove any existing features that overlap any tab region
    df_rec['qualifiers'] = {}
    new_features = []
    for feat in record.features:
        fstart = int(feat.location.start) + 1 + nts_diff
        fend   = int(feat.location.end) - nts_diff
        # check overlap with any tab row
        overlap = False
        for i, row in df_rec.iterrows():
            if not (fend < row['start'] or fstart > row['end']) and feat.type == "CDS":
                overlap = True
                feat.qualifiers['tag'] = [row['tag']]
                feat.qualifiers['product'] = [row['product']]
                feat.qualifiers['mge_element'] = ["yes"]
                new_features.append(feat)
                break
        if not overlap:
            new_features.append(feat)
    record.features = new_features

    # reorder features: keep 'source' first, then CDS (and other) sorted by start
    source_feats = [f for f in record.features if f.type == 'source']
    other_feats = [f for f in record.features if f.type != 'source']
    other_feats.sort(key=lambda f: int(f.location.start))
    record.features = source_feats + other_feats
    return record

def main(biomark, integrases, gbk, output, nts_diff):
    df_biom = load_tab(biomark)
    df_int = pd.read_table(integrases, sep='\t', header=0)
    df_annotation = reformat_table(df_biom, df_int)
    records = list(SeqIO.parse(gbk, 'genbank'))
    annotated = []
    for rec in records:
        annotated.append(annotate_record(rec, df_annotation, nts_diff))
    SeqIO.write(annotated, output, 'genbank')
    print(f"Wrote annotated GenBank to {output}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge Integron_finder biomarkers with GenBank file annotations')
    parser.add_argument('-g', '--gbk', required=True, help='Path to input GenBank file')
    parser.add_argument('-b', '--biomarkers', required=True, help='Path to provirus_genes.tsv file')
    parser.add_argument('-i', '--integrases', required=True, help='Path to provirus.tsv file')
    parser.add_argument('-o', '--output', default='annotated_output.gbk', help='Output GenBank path')
    parser.add_argument('-n', '--nts_diff', default=15, help='Allowed N nucleotide overlapping (default: 15)')
    args = parser.parse_args()

    main(args.biomarkers, args.integrases, args.gbk, args.output, args.nts_diff)
