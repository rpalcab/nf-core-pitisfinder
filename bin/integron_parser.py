#!/usr/bin/env python

# Created by Rosaía Palomino-Cabrera

import pandas as pd
import subprocess
from Bio import SeqIO
from BCBio import GFF
import os
import glob
import sys
import re
import argparse

def get_arguments():

    parser = argparse.ArgumentParser(prog = 'integron_parser.py', description = 'integron_parser.py is part of PITISfinder.')

    input_group = parser.add_argument_group('Input')
    input_group.add_argument('-i', '--integron_file', dest="int_file", required=True, help="Required. Integron file to be filtered and processed", type=os.path.abspath)
    input_group.add_argument('-f', '--fasta_file', dest="fa_file", required=True, help="Required. Fasta file to be processed", type=os.path.abspath)
    input_group.add_argument('-a', '--ann_file', dest="ann_file", required=True, help="Required. Annotation (.gff, .gff3) file to be processed", type=os.path.abspath)
    input_group.add_argument('-r', '--res_file', dest="res_file", required=True, help="Required. Resistance annotation file (.tab) to be processed", type=os.path.abspath)
    input_group.add_argument('-s', '--sample', dest="sample", required=True, help="Required. Sample name", type=os.path.abspath)

    tuning_group = parser.add_argument_group('Tuning parameters')
    tuning_group.add_argument('-n', '--nts_diff', dest="nts_diff", default=10, help="Coordinates differences tolerated when merging annotations. Default: 10", type=int)
    tuning_group.add_argument('-id', '--id_perc', dest="id_perc", default=90, help="Minimum ID percentage required for antomicrobial resistance identification. Default: 90", type=int)
    tuning_group.add_argument('-c', '--q_cov', dest="q_cov", default=80, help="Minimum query coverage required for antomicrobial resistance identification. Default: 80", type=int)

    output_group = parser.add_argument_group('Output')
    output_group.add_argument('-o', '--out_dir', dest='out_dir', required=True, help='Required. Final output folder for reports', type=os.path.abspath)

    arguments = parser.parse_args()

    return arguments

# %%
def read_abr(res_file, contig, d_info, pid, qcov):
    df = pd.read_table(res_file, dtype={'START': 'int64', 'END': 'int64'})
    df.drop(df.loc[(df['%COVERAGE'] < qcov) | (df['%IDENTITY'] < pid)].index, inplace=True)
    df = df[(df.SEQUENCE == contig) & (df.START > d_info['start']) & (df.END < d_info['end'])]
    return df

# %%
def read_gff(ann_file, contig, d_info):
    with open(ann_file) as infile:
        records = list(GFF.parse(infile))
    gff_data = []
    for record in records:
        for feature in record.features:
            gff_data.append([record.id, feature.type, feature.location.start, 
                             feature.location.end, feature.location.strand, feature.qualifiers.get('Name'), feature.qualifiers.get('gene')])

    df0 = pd.DataFrame(gff_data, columns=["seqid", "type", "start", "end", "strand", "name", "gene"])
    df = df0[(df0.type == 'CDS') & (df0.seqid == contig) & (df0.start > d_info['start']) & (df0.end < d_info['end'])]
    return df

def get_annotation(row, df, start, end, annotation, nts_diff):
    df_filtered = df[(abs(df[start] - row['pos_beg']) <= nts_diff) &
                     (abs(df[end] - row['pos_end']) <= nts_diff)]
    if not df_filtered.empty and annotation != None:
        return df_filtered.iloc[0][annotation]
    return None

def merge_annotations(subdf, df_ann, df_abr, nts_diff):
    merged_annotations = []
    for idx, row in subdf.iterrows():
        annot = row.annotation
        if (annot == 'attC') | (annot == 'intI'):
            merged_annotations.append(annot)
            continue
        # Check if present in abricate annotation
        new_annot = get_annotation(row, df_abr, 'START', 'END', 'GENE', nts_diff)
        if pd.isnull(new_annot):
            # If not, check in bakta annotation
            new_annot = get_annotation(row, df_ann, 'start', 'end', 'gene', nts_diff)
        else: 
            merged_annotations.append(new_annot)
            continue
        if pd.isnull(new_annot):
            # Si aún no hay match, usamos la anotación original de subdf
            merged_annotations.append(annot)
        else: 
            merged_annotations.append(new_annot)

    subdf['merged_annotation'] = merged_annotations
    return subdf

# %%
def write_fasta(cds_output_file, d_int):
    with open(cds_output_file, "w") as cds_output_handle:
        for k, v in d_int.items():
            cds_output_handle.write(k + "\n")
            cds_output_handle.write(v + "\n")
    return None

# %%
def extract_fastas(gbk_file, cds_output_file, fna_file, integron):
    # Abrir archivo GBK y los archivos de salida
    with open(gbk_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            flag = 0
            d_int = {}
            d_nucl = {}
            for feature in record.features:
                # Obtener las coordenadas de inicio y final
                start = feature.location.start + 1  # +1 para coordenadas 1-based
                end = feature.location.end
                # Extraer la secuencia del integrón
                if "integron" in feature.type and feature.qualifiers.get('integron_id')[0] == integron:
                    flag = 1
                    integron_seq = feature.extract(record.seq)
                    integron_header = f'>{integron}_{start}_{end}'
                    d_nucl[integron_header] = str(integron_seq)
                # Extraer y guardar las secuencias de CDS
                elif feature.type == "CDS" and flag == 1:
                    # Extraer la secuencia de CDS
                    cds_seq = feature.extract(record.seq)
                    # Obtener ID o locus_tag
                    protein_id = feature.qualifiers.get('protein_id', ['Unknown'])[0]
                    # Escribir en el archivo FASTA
                    header = f">{protein_id}_{start}_{end}"
                    d_int[header] = str(cds_seq)
                
                elif "integron" in feature.type:
                    flag = 0

        # Secuencia nucl. del integrón (fna)
        write_fasta(fna_file, d_nucl)
        # Secuencias CDS del integrón (faa)
        write_fasta(cds_output_file, d_int)
            
    return d_int

# %%
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

    # try:
    #     # extraer faa y fna de los integrones
    #     gbk_file = input_path + f'/{replicon}.gbk'
    #     cds_output_file = input_path + f'/{replicon}_{integron}.faa'
    #     fna_file = input_path + f'/{replicon}_{integron}.fna'
    #     d_int = extract_fastas(gbk_file, cds_output_file, fna_file, integron)
    # except:
    #     print(f'Error en el parseo de GBK. Comprueba la integridad de {gbk_file}')
    #     return None

    # df_abr = pd.DataFrame(columns = ['pos_beg', 'pos_end', 'abr_ann'])
    # df_prokka = pd.DataFrame(columns = ['pos_beg', 'pos_end', 'prokka_ann'])
    # # Si hay CDS en el integrón, anotarlos con prokka y abricate
    # if len(d_int.keys()) > 0:
    #     prokka_dir, abr_out = annotate_cds(cds_output_file, input_path, replicon, integron)
    #     df_prokka = prokka_parse(prokka_dir)
    #     df_abr = abr_parse(abr_out)
    # # mergear con las anotaciones, priorizar abricate y reorientar df si necesario
    # subdf['pos_beg'] = subdf['pos_beg'].astype('int64')
    # subdf['pos_end'] = subdf['pos_end'].astype('int64')
    # mid_df = pd.merge(subdf, df_prokka, on=['pos_beg', 'pos_end'], how='outer')
    # final_df = pd.merge(mid_df, df_abr, on=['pos_beg', 'pos_end'], how='outer')
    # final_df['ann'] = final_df['abr_ann'].fillna(final_df['prokka_ann'])
    # final_df.drop_duplicates(inplace=True)

    
    # # group cassettes
    # cassettes = []
    # current_cassette = []

    # for _, row in final_df.iterrows():
    #     if row['type_elt'] == 'attC':
    #         cassettes.append(current_cassette)
    #         current_cassette = []
    #     elif row['type_elt'] == 'protein' and row['annotation'] != 'intI':
    #         if pd.isnull(row['ann']):
    #             ann = "NA;hypothetical protein"
    #         else:
    #             ann = row['ann']
    #         current_cassette.append(ann)

    # # Por si acaso no se reconoce el último attC
    # cassettes.append(current_cassette)
    # # y eliminamos cassettes vacíos (?)
    # cassettes = [i for i in cassettes if i != []]
    # genes = '_'.join([i[0].split(';')[0] for i in cassettes])
    # genes = re.sub(r'[^a-zA-Z0-9\-\_]', '', genes)


    # name = f'int_{genes}_{sample}_{count}'
    # # for i in cassettes: info.append(i)
    # for i in cassettes: info.append(','.join(i))
    # info.extend([""] * (20-len(info)))

    # # Guardamos secuencia nucleotídica
    # cp_cmd = ['cp', fna_file, f'{original_path}/11_integrons/{name}.fasta']
    # subprocess.run(cp_cmd)

# %%
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

    # Create empty summary df
    summary_df = pd.DataFrame(columns=['Sample', 'Contig', 'Name', 'Size', 'Start', 'End', 'Integrase', 'Cassette 1',
                                        'Cassette 2', 'Cassette 3', 'Cassette 4', 'Cassette 5', 'Cassette 6', 'Cassette 7', 
                                        'Cassette 8', 'Cassette 9', 'Cassette 10', 'Cassette 11', 'Cassette 12'])

    # Divide into integrons and chromosomes
    df_grouped = df_integrons.groupby(['ID_replicon', 'ID_integron'])

    for count, ((contig, integron), subdf) in enumerate(df_grouped):    # keep count of integron number
        d_info, subdf = basic_info(subdf)
        df_ann = read_gff(ann_file, contig, d_info)
        df_abr = read_abr(res_file, contig, d_info, args.id_perc, args.q_cov)
        df_ann.to_csv(f'{report_out}_gff')
        unified_df = merge_annotations(subdf, df_ann, df_abr, args.nts_diff)
        # summary_df.loc[len(summary_df)] = info

    summary_df.to_csv(report_out, index=False)


    # Integron files

    # original_path = os.path.abspath(sys.argv[1])
    # summary_df = pd.DataFrame(columns=['Sample', 'Pl/Chr', 'Name', 'Size', 'Inicio', 'Final', 'Tipo', 'Integrasa', 'Cassette 1',
    #                                     'Cassette 2', 'Cassette 3', 'Cassette 4', 'Cassette 5', 'Cassette 6',
    #                                     'Cassette 7', 'Cassette 8', 'Cassette 9', 'Cassette 10', 'Cassette 11', 'Cassette 12'])
    # for integron_file in glob.glob(f'{original_path}/11_integrons/*/Results_Integron_Finder_*/*.integrons'):
    #     input_path = os.path.dirname(os.path.abspath(integron_file))
    #     sample = input_path.split('/')[-2]
    #     print(sample)
    #     try:
    #         df_integron = pd.read_table(integron_file, comment='#')
    #     except:
    #         print(f'No integrons in {sample}')
    #         continue

    #     # Divide into integrons and chromosomes (cada cromosoma resetea el número de integrón, así que puede haber varios integron_01 por muestra)
    #     grouped = df_integron.groupby(['ID_replicon', 'ID_integron'])
    #     subdfs = {}

    #     for count, ((replicon, integron), group) in enumerate(grouped):
    #         key = f"{replicon}_{integron}"
    #         subdfs[key] = group
    #         info = extract_info(sample, subdfs[key], replicon, integron, input_path, original_path, count)
    #         if info:
    #             summary_df.loc[len(summary_df)] = info

    # summary_df.to_csv(f'{original_path}/11_integrons/integron_summary.csv', index=False)
    # print(f'Integron summary in {original_path}/11_integrons/integron_summary.csv')
