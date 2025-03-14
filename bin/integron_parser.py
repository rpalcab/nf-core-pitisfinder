#!/usr/bin/env python

# Created by Rosaía Palomino-Cabrera

import pandas as pd
from Bio import SeqIO
from BCBio import GFF
import os
import re
import argparse

def get_arguments():

    parser = argparse.ArgumentParser(prog = 'integron_parser.py', description = 'integron_parser.py is part of PITISfinder.')

    input_group = parser.add_argument_group('Input')
    input_group.add_argument('-i', '--integron_file', dest="int_file", type=os.path.abspath,
                             required=True, help="Required. Integron file to be filtered and processed")
    input_group.add_argument('-f', '--fasta_file', dest="fa_file", required=True, type=os.path.abspath, 
                             help="Required. Fasta file to be processed")
    input_group.add_argument('-a', '--ann_file', dest="ann_file", required=True, type=os.path.abspath, 
                             help="Required. Annotation (.gff, .gff3) file to be processed")
    input_group.add_argument('-r', '--res_file', dest="res_file", required=True, type=os.path.abspath, 
                             help="Required. Resistance annotation file (.tab) to be processed")
    input_group.add_argument('-s', '--sample', dest="sample", required=True, help="Required. Sample name", type=str)

    tuning_group = parser.add_argument_group('Tuning parameters')
    tuning_group.add_argument('--nts_diff', dest="nts_diff", default=10, type=int, 
                              help="Coordinates differences tolerated when merging annotations. Default: 10")
    tuning_group.add_argument('--id_perc', dest="id_perc", default=90, type=int, 
                              help="Minimum ID percentage required for antomicrobial resistance identification. Default: 90")
    tuning_group.add_argument('--q_cov', dest="q_cov", default=80, type=int, 
                              help="Minimum query coverage required for antomicrobial resistance identification. Default: 80")
    tuning_group.add_argument('--max_cas', dest="max_cassettes", default=12, type=int, 
                              help="Maximum cassettes in integron. Default: 12")
    
    output_group = parser.add_argument_group('Output')
    output_group.add_argument('-o', '--out_dir', dest='out_dir', required=True, type=os.path.abspath,
                              help='Required. Final output folder for reports')

    arguments = parser.parse_args()

    return arguments

def basic_info(subdf):
    d_info = {}
    # Integrase info
    integrase_row = subdf[subdf['annotation'] == 'intI'].iloc[0]
    d_info['integrase_model'] = integrase_row['model']
    integrase_strand = integrase_row['strand']
    if integrase_strand > 0:        # Redirect dataframe (int is always at the beginning)
        subdf = subdf[::-1]
    # Integron coords and size
    d_info['start'] = min(subdf['pos_beg'])
    d_info['end'] = max(subdf['pos_end'])
    d_info['size'] = d_info['end'] - d_info['start']
    return d_info, subdf

def read_gff(ann_file, contig, d_info):
    with open(ann_file) as infile:                                                              # Open gff file and parse
        records = list(GFF.parse(infile))
    gff_data = []
    for record in records:
        for feature in record.features:
            gff_data.append([record.id, feature.type, feature.location.start, 
                             feature.location.end, feature.location.strand, 
                             feature.qualifiers.get('Name', [None])[0], feature.qualifiers.get('gene', [None])[0]])

    df0 = pd.DataFrame(gff_data, columns=["seqid", "type", "start", "end", "strand", "name", "gene"])
    df = df0[(df0.type == 'CDS') & (df0.seqid == contig) & 
             (df0.start > d_info['start']) & (df0.end < d_info['end'])]   # Keep genes inside integron
    return df

def read_abr(res_file, contig, d_info, pid, qcov):
    df = pd.read_table(res_file, dtype={'START': 'int64', 'END': 'int64'})
    df.drop(df.loc[(df['%COVERAGE'] < qcov) | (df['%IDENTITY'] < pid)].index, inplace=True)     # Remove low quality hits
    df = df[(df.SEQUENCE == contig) & (df.START > d_info['start']) & (df.END < d_info['end'])]  # Keep genes inside integron
    return df

def read_annotation(row, df, start, end, annotation, nts_diff):
    df_filtered = df[(abs(df[start] - row['pos_beg']) <= nts_diff) &    # Takes annotation in region
                     (abs(df[end] - row['pos_end']) <= nts_diff)]
    return df_filtered[annotation].iloc[0] if not df_filtered.empty else None   # Returns annotation if it exists

def update_annotation(row, df_ann, df_abr, nts_diff):
        if row.annotation in {'attC', 'intI'}:
            return row.annotation
        
        new_annot = read_annotation(row, df_abr, 'START', 'END', 'GENE', nts_diff) or \
                    read_annotation(row, df_ann, 'start', 'end', 'gene', nts_diff)
        
        return new_annot if new_annot is not None else row.annotation

def merge_annotations(subdf, df_ann, df_abr, nts_diff):
    subdf = subdf.copy()
    subdf['merged_annotation'] = subdf.apply(lambda row:update_annotation(row, df_ann, df_abr, nts_diff), axis=1)
    return subdf

def get_cassettes(unified_df, count, sample):
    l_cassettes = unified_df.loc[(unified_df['type_elt'] == 'protein') & 
                                 (unified_df['annotation'] != 'intI')]['merged_annotation'].to_list()
    cassettes = '_'.join([re.sub(r'[^a-zA-Z0-9]', '', i) for i in l_cassettes])
    name = f'int_{cassettes}_{sample}_{count}'
    return l_cassettes, name

def merge_info(summary_df, sample, contig, integron_name, d_info, cassettes, max_cassettes):
    cassettes = (cassettes[:max_cassettes] + [None] * 
                (max_cassettes - len(cassettes[:max_cassettes])))                             # Truncate cassettes if too many, add None if too few
    row = [sample, contig, integron_name, d_info['size'], d_info['start'], 
           d_info['end'], d_info['integrase_model']] + cassettes    # Create row to append
    summary_df.loc[len(summary_df)] = row
    return summary_df

def write_fasta(cds_output_file, d_int):
    with open(cds_output_file, "w") as cds_output_handle:
        for k, v in d_int.items():
            cds_output_handle.write(f'>{k}' + "\n")
            cds_output_handle.write(v + "\n")
    return None

def save_sequence(unified_df, fa_file):
    cont = unified_df['Contig']
    start = unified_df['Start'] - 1
    end = unified_df['End'] - 1
    name = unified_df['Name']

    for seq_record in SeqIO.parse(fa_file, "fasta"):
        if seq_record.id == cont:
            d = {name: str(seq_record.seq[start:end])}
            write_fasta(f'{name}.fasta', d)
            break
    return None



if __name__ == "__main__":

    # Parameters
    args = get_arguments()

    # Input
    int_file = args.int_file
    fa_file = args.fa_file
    ann_file = args.ann_file
    res_file = args.res_file
    sample = args.sample

    # Output
    out_dir = args.out_dir
    report_out = os.path.join(out_dir, "integrons_summary.tsv")

    # Try open integron file, quit if fail
    try:
        df_integrons = pd.read_table(int_file, comment='#')
    except:
        print(f'No integrons in {sample}')
        exit()

    # Keep only complete integrons, quit if none
    df_integrons = df_integrons[df_integrons.type == 'complete']
    if len(df_integrons) == 0:
        print(f'No complete integrons in {sample}')
        exit()

    # Empty dataframe with {max_cassettes} cassettes
    summary_df = pd.DataFrame(columns=['Sample', 'Contig', 'Name', 'Size', 'Start', 'End', 'Integrase'] +   
                                      [f'Cassette {i+1}' for i in range(args.max_cassettes)])
    
    # Divide into integrons and chromosomes
    df_grouped = df_integrons.groupby(['ID_replicon', 'ID_integron'])

    for count, ((contig, integron), subdf) in enumerate(df_grouped):    # keep count of integron number
        d_info, subdf = basic_info(subdf)                               # Read coordinates
        df_ann = read_gff(ann_file, contig, d_info)                     # Read bakta annotation
        df_abr = read_abr(res_file, contig, d_info, args.id_perc, args.q_cov)   # Read abricate annotation
        unified_df = merge_annotations(subdf, df_ann, df_abr, args.nts_diff)    # Update annotations: abr > bakta > i_f
        cassettes, integron_name = get_cassettes(unified_df, count, sample)                     # get cassette, get name
        summary_df = merge_info(summary_df, sample, contig, integron_name, d_info, cassettes, args.max_cassettes)
        save_sequence(summary_df.iloc[-1], fa_file)                                    # Save integron sequence (read fasta, get coordinates, save sequence)

    summary_df.to_csv(report_out, index=False)
