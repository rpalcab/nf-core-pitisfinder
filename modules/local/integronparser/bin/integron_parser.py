# %%
import pandas as pd
import subprocess
from Bio import SeqIO
from BCBio import GFF
import os
import glob
import sys
import re

# %%
def abr_parse(abr_out):
    abr_raw = pd.read_table(abr_out)
    if len(abr_raw) > 0:
        abr_raw[['pos_beg', 'pos_end']] = pd.DataFrame(abr_raw['SEQUENCE'].str.split('_').str[-2:].tolist(), index=abr_raw.index)
        abr_raw.drop(abr_raw.loc[abr_raw['%IDENTITY'] < 90].index, inplace=True)
        abr_raw.drop(abr_raw.loc[abr_raw['%COVERAGE'] < 80].index, inplace=True)
        df_abr = abr_raw[['pos_beg', 'pos_end', 'GENE']]
        df_abr.columns = ['pos_beg', 'pos_end', 'abr_ann']
        df_abr['pos_beg'] = df_abr['pos_beg'].astype('int64')
        df_abr['pos_end'] = df_abr['pos_end'].astype('int64')
        return df_abr
    return pd.DataFrame(columns=['pos_beg', 'pos_end', 'abr_ann'])

# %%
def prokka_parse(prokka_dir):
    df_prokka = pd.DataFrame(columns=['pos_beg', 'pos_end', 'prokka_ann'])
    for prokka_file in glob.glob(f'{prokka_dir}/*.gff'):
        in_handle = open(prokka_file)
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                start, end = rec.id.split('_')[-2:]
                gene = feature.qualifiers.get('gene',['NA'])[0]
                product = feature.qualifiers.get('product')[0]
                df_prokka.loc[-1] = [start, end, f'{gene};{product}']
                df_prokka.index = df_prokka.index + 1
                df_prokka = df_prokka.sort_index()
        in_handle.close()
    df_prokka['pos_beg'] = df_prokka['pos_beg'].astype('int64')
    df_prokka['pos_end'] = df_prokka['pos_end'].astype('int64')
    return df_prokka

# %%
def write_fasta(cds_output_file, d_int):
    with open(cds_output_file, "w") as cds_output_handle:
        for k, v in d_int.items():
            cds_output_handle.write(k + "\n")
            cds_output_handle.write(v + "\n")
    return None

# %%
def annotate_cds(cds_output_file, input_path, replicon, integron):
    # anotación con prokka
    prokka_cmd = ['prokka', '--quiet', '--force', cds_output_file, '--outdir', input_path + f'/prokka_{replicon}_{integron}']
    subprocess.run(prokka_cmd)
    # anotación con abricate
    abr_cmd1 = ['abricate', cds_output_file]
    abr_out = open(input_path + f'/abricate_{replicon}_{integron}.out', 'w')
    subprocess.run(abr_cmd1, stdout=abr_out)
    return input_path + f'/prokka_{replicon}_{integron}', input_path + f'/abricate_{replicon}_{integron}.out'

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
def extract_info(sample, subdf, replicon, integron, input_path, original_path, count):
    # integrase
    try:
        integrase_row = subdf[subdf['annotation'] == 'intI'].iloc[0]
        integrase_model = integrase_row['model']
        integrase_strand = integrase_row['strand']
    except:
        integrase_model = ''
        integrase_strand = -1

    # attC sites
    try:
        attc_models = subdf[subdf['type_elt'] == 'attC']['model'].tolist()
    except:
        attc_models = ''

    try:
        # extraer faa y fna de los integrones
        gbk_file = input_path + f'/{replicon}.gbk'
        cds_output_file = input_path + f'/{replicon}_{integron}.faa'
        fna_file = input_path + f'/{replicon}_{integron}.fna'
        d_int = extract_fastas(gbk_file, cds_output_file, fna_file, integron)
    except:
        print(f'Error en el parseo de GBK. Comprueba la integridad de {gbk_file}')
        return None

    df_abr = pd.DataFrame(columns = ['pos_beg', 'pos_end', 'abr_ann'])
    df_prokka = pd.DataFrame(columns = ['pos_beg', 'pos_end', 'prokka_ann'])
    # Si hay CDS en el integrón, anotarlos con prokka y abricate
    if len(d_int.keys()) > 0:
        prokka_dir, abr_out = annotate_cds(cds_output_file, input_path, replicon, integron)
        df_prokka = prokka_parse(prokka_dir)
        df_abr = abr_parse(abr_out)
    # mergear con las anotaciones, priorizar abricate y reorientar df si necesario
    subdf['pos_beg'] = subdf['pos_beg'].astype('int64')
    subdf['pos_end'] = subdf['pos_end'].astype('int64')
    mid_df = pd.merge(subdf, df_prokka, on=['pos_beg', 'pos_end'], how='outer')
    final_df = pd.merge(mid_df, df_abr, on=['pos_beg', 'pos_end'], how='outer')
    final_df['ann'] = final_df['abr_ann'].fillna(final_df['prokka_ann'])
    final_df.drop_duplicates(inplace=True)
    if integrase_strand > 0:
        final_df = final_df[::-1]
    
    # group cassettes
    cassettes = []
    current_cassette = []

    for _, row in final_df.iterrows():
        if row['type_elt'] == 'attC':
            cassettes.append(current_cassette)
            current_cassette = []
        elif row['type_elt'] == 'protein' and row['annotation'] != 'intI':
            if pd.isnull(row['ann']):
                ann = "NA;hypothetical protein"
            else:
                ann = row['ann']
            current_cassette.append(ann)

    # Por si acaso no se reconoce el último attC
    cassettes.append(current_cassette)
    # y eliminamos cassettes vacíos (?)
    cassettes = [i for i in cassettes if i != []]
    genes = '_'.join([i[0].split(';')[0] for i in cassettes])
    genes = re.sub(r'[^a-zA-Z0-9\-\_]', '', genes)

    # lista con datos para Output final
    contig = final_df['ID_replicon'].values[0]
    start = min(final_df['pos_beg'])
    end = max(final_df['pos_end'])
    size = end - start
    type = final_df['type'].values[0]
    if type == 'In0':
        name = f'In0_{sample}_{count}'
    else:
    # name = f'{integron}_{replicon}_{genes}_{sample}'
    # Int_AMR_AMR_..._Idunico
        name = f'int_{genes}_{sample}_{count}'

    info = [sample, contig, name, size, start, end, type, integrase_model]
    # for i in cassettes: info.append(i)
    for i in cassettes: info.append(','.join(i))
    info.extend([""] * (20-len(info)))

    # Guardamos secuencia nucleotídica
    cp_cmd = ['cp', fna_file, f'{original_path}/11_integrons/{name}.fasta']
    subprocess.run(cp_cmd)

    return info

# %%
if __name__ == "__main__":
    # Integron files
    original_path = os.path.abspath(sys.argv[1])
    summary_df = pd.DataFrame(columns=['Sample', 'Pl/Chr', 'Name', 'Size', 'Inicio', 'Final', 'Tipo', 'Integrasa', 'Cassette 1',
                                        'Cassette 2', 'Cassette 3', 'Cassette 4', 'Cassette 5', 'Cassette 6',
                                        'Cassette 7', 'Cassette 8', 'Cassette 9', 'Cassette 10', 'Cassette 11', 'Cassette 12'])
    for integron_file in glob.glob(f'{original_path}/11_integrons/*/Results_Integron_Finder_*/*.integrons'):
        input_path = os.path.dirname(os.path.abspath(integron_file))
        sample = input_path.split('/')[-2]
        print(sample)
        try:
            df_integron = pd.read_table(integron_file, comment='#')
        except:
            print(f'No integrons in {sample}')
            continue

        # Divide into integrons and chromosomes (cada cromosoma resetea el número de integrón, así que puede haber varios integron_01 por muestra)
        grouped = df_integron.groupby(['ID_replicon', 'ID_integron'])
        subdfs = {}

        for count, ((replicon, integron), group) in enumerate(grouped):
            key = f"{replicon}_{integron}"
            subdfs[key] = group
            info = extract_info(sample, subdfs[key], replicon, integron, input_path, original_path, count)
            if info:
                summary_df.loc[len(summary_df)] = info

    summary_df.to_csv(f'{original_path}/11_integrons/integron_summary.csv', index=False)
    print(f'Integron summary in {original_path}/11_integrons/integron_summary.csv')
