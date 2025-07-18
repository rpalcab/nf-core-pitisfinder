
#!/usr/bin/env python3

# %%
import argparse
import pandas as pd
import glob
import plotly.express as px
import plotly.graph_objects as go
from jinja2 import Template
import numpy as np
import math
import shutil
from pathlib import Path

# %%
template_str = """
{% macro render_filter_cell(col, col_type, filter_options) %}
    {% if col_type == 'numeric' %}
        <th>
            <div class="min-max-wrapper">
                <input type="number" class="min filter-box" placeholder="Min">
                <input type="number" class="max filter-box" placeholder="Max">
            </div>
        </th>
    {% elif col_type == 'text' %}
        <th><input type="text" class="filter-box" placeholder="Search"></th>
    {% elif col_type == 'categorical' %}
        <th>
            <select class="filter-box">
                <option value="">All</option>
                {% for val in filter_options.get(col, []) %}
                    <option value="{{ val }}">{{ val }}</option>
                {% endfor %}
            </select>
        </th>
    {% else %}
        <th></th>
    {% endif %}
{% endmacro %}

{% macro render_table(df, zipped_columns, filter_options, table_id) %}
<div style="overflow-x: auto;">
    <table id="{{ table_id }}" class="display nowrap">
        <thead>
            <tr>
                {% for col, _ in zipped_columns %}
                    <th>{{ col }}</th>
                {% endfor %}
            </tr>
            <tr class="filter-row">
                {% for col, col_type in zipped_columns %}
                    {{ render_filter_cell(col, col_type, filter_options) }}
                {% endfor %}
            </tr>
        </thead>
        <tbody>
            {% for row in df.itertuples(index=False) %}
                <tr>
                    {% for cell in row %}
                        <td>{{ cell | safe }}</td>
                    {% endfor %}
                </tr>
            {% endfor %}
        </tbody>
    </table>
</div>
{% endmacro %}

{% macro render_section(title, table_html, table_id, figs=[]) %}
<details class="section" open>
    <summary>
        <span class="chevron"></span>
        <h2>{{ title }}</h2>
    </summary>
    {{ table_html }}
    <div id="sample-filter-container-{{ table_id }}" class="sample-filter-container" style="margin-top: 10px;"></div>
    {% if figs %}
    <div class="plot-row" id="plot_{{ table_id }}">
        {% for fig in figs %}
            <div>{{ fig | safe }}</div>
        {% endfor %}
    </div>
    {% endif %}
</details>
{% endmacro %}

{% macro render_empty(title, mge) %}
<details class="section" open>
    <summary>
        <span class="chevron"></span>
        <h2>{{ title }}</h2>
    </summary>
    <p>No {{ mge }} found.</p>
</details>
{% endmacro %}

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>MGE full report</title>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/buttons/2.4.1/css/buttons.dataTables.min.css">
    <link href="https://cdn.jsdelivr.net/npm/select2@4.1.0-rc.0/dist/css/select2.min.css" rel="stylesheet" />
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css" />
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .section { margin-bottom: 80px; }
        h2 { border-bottom: 1px solid #ccc; padding-bottom: 6px; }
        .dataTables_wrapper { width: 100%; overflow-x: auto; }
        table.dataTable { table-layout: auto; width: 100%; word-break: break-word; border-collapse: collapse; }
        table.dataTable tbody tr:nth-child(odd) { background-color: #f9f9f9; }
        table.dataTable tbody tr:nth-child(even) { background-color: #ffffff; }
        table.dataTable th, table.dataTable td { border: none; padding: 6px 10px; text-align: center; word-break: break-word; }
        thead tr.filter-row th { padding: 4px; }
        .filter-box { width: 100%; box-sizing: border-box; }
        .min-max-wrapper { display: flex; gap: 4px; }
        .min-max-wrapper input { width: 48%; }
        .dt-left-toolbar select { min-width: 100px; margin-right: 0.5rem; }
        input[type='number']::-webkit-outer-spin-button, input[type='number']::-webkit-inner-spin-button { -webkit-appearance: none; margin: 0; }
        input[type='number'] { -moz-appearance: textfield; }
        .plot-row { display: flex; flex-wrap: wrap; gap: 40px; }
        .plot-row iframe { flex: 1 1 48%; min-width: 300px; margin: 0; padding: 0; }
        details.section { margin-bottom: 10px; border: 1px solid #ccc; border-radius: 6px; padding: 10px 14px; background-color: #ffffff; }
        details > summary { cursor: pointer; display: flex; align-items: center; list-style: none; outline: none; }
        details > summary:hover { background-color: #ffffff; }
        details > summary h2 { display: inline; margin: 0; padding-left: 8px; font-size: 1.2em; }
        details > summary::-webkit-details-marker { display: none; }
        .chevron { display: inline-block; width: 12px; height: 12px; margin-right: 8px; border-right: 2px solid #333; border-bottom: 2px solid #333; transform: rotate(45deg); transition: transform 0.2s ease; margin-top: 2px; }
        details[open] > summary .chevron { transform: rotate(135deg); margin-top: 4px; }

        // REMOVE
        .dt-left-toolbar { background-color: #e8f4ff; border: 1px solid #99c; }
        .dt-right-toolbar { background-color: #ffe8e8; border: 1px solid #c99; }

    </style>
</head>
<body>

{{ render_section("1. General Report", render_table(df_general, zipped_columns_general, filter_options_general, 'table_general'), 'table_general', [d_figs.get('general')]) }}

{% if df_pl is not none %}
    {{ render_section("2. Plasmids", render_table(df_pl, zipped_columns_pl, filter_options_pl, 'table_pl'), 'table_pl', d_figs.get('plasmid', [])) }}
{% else %}
    {{ render_empty('2. Plasmids', 'plasmids') }}
{% endif %}

{% if df_int is not none %}
    {{ render_section("3. Integrons", render_table(df_int, zipped_columns_int, filter_options_int, 'table_int'), 'table_int', d_figs.get('integron', [])) }}
{% else %}
    {{ render_empty('3. Integrons', 'integrons') }}
{% endif %}

{% if df_ph is not none %}
    {{ render_section("4. Prophages", render_table(df_ph, zipped_columns_ph, filter_options_ph, 'table_ph'), 'table_ph', d_figs.get('prophage', [])) }}
{% else %}
    {{ render_empty('4. Prophages', 'prophages') }}
{% endif %}

{% if df_is is not none %}
    {{ render_section("5. IS", render_table(df_is, zipped_columns_is, filter_options_is, 'table_is'), 'table_is', d_figs.get('IS', [])) }}
{% else %}
    {{ render_empty('5. IS', 'IS') }}
{% endif %}

<script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
<script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.1/js/dataTables.buttons.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.1/js/buttons.colVis.min.js"></script>

<script>
    function initTableFilters(tableId, colTypes) {
        // Build an array of indexes for all numeric columns
        const numericCols = colTypes
            .map((t, i) => t === 'numeric' ? i : null)
            .filter(i => i !== null);

        const table = $('#' + tableId).DataTable({
            orderCellsTop: true,
            fixedHeader: true,
            pageLength: 10,
            autoWidth: false,
            scrollX: true,
            dom:
                "<'row align-items-center mb-2'<'col-sm-12 d-flex align-items-center dt-left-toolbar'<'filter-sample me-3'>Bf>>" +
                "<'row'<'col-sm-12'tr>>" +
                "<'row'<'col-sm-5'l><'col-sm-2'i><'col-sm-5'p>>",
            buttons: [
                {
                    extend: 'colvis',
                    text: 'Select columns',
                    className: 'me-2'
                }
            ],
            initComplete: function () {
                const api = this.api();

                // 1) Sample‑column dropdown
                const sampleColIndex = 0;
                const uniqueSamples = api.column(sampleColIndex).data().unique().sort();
                const select = $(
                    '<select class="form-select form-select-sm ms-2" style="width: auto;">' +
                    '<option value="">All Samples</option>' +
                    '</select>'
                );
                const container = api.table().container();
                uniqueSamples.each(d => {
                    if (d) select.append(`<option value="${d}">${d}</option>`);
                });
                $(api.table().container()).find('.dt-left-toolbar')
                    .append($('<label class="me-2 mb-0">Filter Sample:</label>'))
                    .append(select);
                select.on('change', () => {
                    const v = select.val();
                    api.column(sampleColIndex)
                       .search(v ? '^' + v + '$' : '', true, false)
                       .draw();
                });

                // 2) Numeric‑range filter
                $.fn.dataTable.ext.search.push((settings, data) => {
                    if (settings.nTable.id !== tableId) return true;

                    for (const i of numericCols) {
                        const rawMin = $( container )
                            .find('tr.filter-row th').eq(i)
                            .find('input.min').val();
                        const rawMax = $( container )
                            .find('tr.filter-row th').eq(i)
                            .find('input.max').val();
                        const min    = rawMin === '' ? -Infinity : parseFloat(rawMin);
                        const max    = rawMax === '' ?  Infinity : parseFloat(rawMax);
                        const cell   = parseFloat(data[i]);

                        // Only reject if a bound is specified and the cell is outside it
                        if ((rawMin !== '' && cell < min) ||
                            (rawMax !== '' && cell > max)) {
                            return false;
                        }
                    }
                    return true;
                });

                // 3) Per‑column text/categorical inputs
                api.columns().every(function (i) {
                    const header = $(this.header());
                    const input  = header.closest('thead')
                                         .find('tr.filter-row th')
                                         .eq(i)
                                         .find('input, select');
                    if (!input.length) return;

                    input.on('input change', () => {
                        const v = input.val() || '';
                        if (colTypes[i] === 'numeric') {
                            table.draw();           // redraw to trigger ext.search
                        } else {
                            api.column(i).search(v, false, false).draw();
                        }
                    });
                });
            }
        });
    }

    $(document).ready(function () {
        initTableFilters('table_general', {{ column_types_general | tojson }});
        {% if df_pl is not none %}initTableFilters('table_pl', {{ column_types_pl | tojson }});{% endif %}
        {% if df_int is not none %}initTableFilters('table_int', {{ column_types_int | tojson }});{% endif %}
        {% if df_ph is not none %}initTableFilters('table_ph', {{ column_types_ph | tojson }});{% endif %}
        {% if df_is is not none %}initTableFilters('table_is', {{ column_types_is | tojson }});{% endif %}
    });
</script>


</body>
</html>
"""

# %%
def infer_column_types(df):
    col_types = []
    for col in df.columns:
        series = df[col].replace("", np.nan)
        if pd.api.types.is_numeric_dtype(series) and not all(series.isna()):
            col_types.append("numeric")
        elif df[col].nunique() < 5 and col not in ['AMR', 'VF', 'DF']:
            col_types.append("categorical")
        elif col == 'Preview':
            col_types.append("preview")
        else:
            col_types.append("text")
    return col_types

# %%
def inject_preview_icons(df):
    icon_html = '<i class="fas fa-image"></i>'

    df = df.copy()
    df['Preview'] = df['Preview'].apply(
        lambda url: f'<a href="{url}" target="_blank" rel="noopener">{icon_html}</a>'
    )
    return df

# %%
def prepare_table_data(df):
    if df.empty:
        return None, None, None
    col_types = infer_column_types(df)
    zipped = list(zip(df.columns, col_types))
    filters = {
        col: sorted(df[col].dropna().astype(str).unique())
        for col, t in zip(df.columns, col_types)
        if t == "categorical"
    }
    return col_types, zipped, filters

# %%
def render_final_report(df_general, df_pl, df_int, df_ph, df_is, d_figs, input_dir, output_file):
    outpath = Path(output_file).parent / "images"
    outpath.mkdir(parents=True, exist_ok=True)

    # Copy images to output directory
    df_general['tmp_img'] = df_general.apply(lambda row: row['Name'] if row['MGE'] != 'plasmid' else f'{row['Name'].split(':')[0]}_{row['Sample']}', axis=1)
    df_general.apply(lambda row: shutil.copy2(f'{input_dir}/{row['Sample']}/{row['MGE']}s/{row['tmp_img']}.png', outpath), axis=1)
    df_general['Preview'] = df_general.apply(lambda row: f'{outpath}/{row['tmp_img']}.png', axis=1)
    df_general.drop(columns=['tmp_img'], inplace=True)
    df_general = inject_preview_icons(df_general)
    # Copy to mge dataframes
    basic_df = df_general[['Sample', 'Contig', 'Name', 'MGE', 'Preview']]
    df_pl = df_pl.merge(basic_df[basic_df['MGE'] == 'plasmid'][['Sample', 'Contig', 'Preview']], on=['Sample', 'Contig'], how='left')
    df_int = df_int.merge(basic_df[basic_df['MGE'] == 'integron'][['Sample', 'Contig', 'Name', 'Preview']], on=['Sample', 'Contig', 'Name'], how='left').replace(np.nan, '')
    df_ph = df_ph.merge(basic_df[basic_df['MGE'] == 'prophage'][['Sample', 'Contig', 'Name', 'Preview']], on=['Sample', 'Contig', 'Name'], how='left')

    column_types_general, zipped_columns_general, filter_options_general = prepare_table_data(df_general)
    column_types_pl, zipped_columns_pl, filter_options_pl = prepare_table_data(df_pl)
    column_types_int, zipped_columns_int, filter_options_int = prepare_table_data(df_int)
    column_types_ph, zipped_columns_ph, filter_options_ph = prepare_table_data(df_ph)
    column_types_is, zipped_columns_is, filter_options_is = prepare_table_data(df_is)

    rendered = Template(template_str).render(
        df_general=df_general,
        df_pl=df_pl if not df_pl.empty else None,
        df_int=df_int if not df_int.empty else None,
        df_ph=df_ph if not df_ph.empty else None,
        df_is=df_is if not df_is.empty else None,
        d_figs=d_figs,
        column_types_general=column_types_general,
        zipped_columns_general=zipped_columns_general,
        filter_options_general=filter_options_general,
        column_types_pl=column_types_pl,
        zipped_columns_pl=zipped_columns_pl,
        filter_options_pl=filter_options_pl,
        column_types_int=column_types_int,
        zipped_columns_int=zipped_columns_int,
        filter_options_int=filter_options_int,
        column_types_ph=column_types_ph,
        zipped_columns_ph=zipped_columns_ph,
        filter_options_ph=filter_options_ph,
        column_types_is=column_types_is,
        zipped_columns_is=zipped_columns_is,
        filter_options_is=filter_options_is,
    )
    
    with open(output_file, "w") as f:
        f.write(rendered)

# %%
def load_concat(reports):
    if len(reports) == 0:
        return pd.DataFrame()
    l_dfs = [pd.read_csv(report, sep='\t', header=0) for report in reports]
    df = pd.concat(l_dfs)
    df.fillna('', inplace=True)
    if "Start" in df.columns:
        df["Start"] = df["Start"].apply(lambda x: int(x) if x != "" else "")
        df["End"] = df["End"].apply(lambda x: int(x) if x != "" else "")

    return df

# %%
def plot_general(df_general):
    df_counts = df_general.groupby(['Sample', 'MGE']).size().reset_index(name='Count')

    fig = px.bar(
        df_counts,
        x='Sample',
        y='Count',
        color='MGE',
        text='Count',
        labels={'Count': 'Number of MGEs'},
        category_orders={'Sample': sorted(df_general['Sample'].unique())}  # consistent sample order
    )

    fig.update_layout(
        barmode='stack',
        xaxis_title='Sample',
        yaxis_title='MGE count',
        title={
            'text': 'Count of MGEs per sample',
            'y':1,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        width=1200,
        height=400,
    )
    fig.update_traces(textposition='auto', insidetextanchor='middle', cliponaxis=True)

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

# %%
def plot_size(mge, df_sub):
    jitter_strength = 0.1
    x_center = 0
    x_jittered = np.random.normal(loc=x_center, scale=jitter_strength, size=len(df_sub))

    fig = go.Figure()

    fig.add_trace(go.Box(
        y=df_sub['Length_log'],
        x=[x_center] * len(df_sub),
        marker_color='#2b9eb3',
        name=mge,
        boxpoints=False,
    ))

    fig.add_trace(go.Scatter(
        y=df_sub['Length_log'],
        x=x_jittered,
        mode='markers',
        marker=dict(size=6, color='black', opacity=0.5),
        customdata=np.stack((df_sub['Name'], df_sub['Length']), axis=-1),
        hovertemplate='Name: %{customdata[0]}<br>Length: %{customdata[1]} bp<extra></extra>',
        name=''
    ))

    fig.update_layout(
        title={
            'text': f'Length distribution in {mge}s',
            'y':1,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        yaxis_title='Length (bp)' + (' [log scale]' if mge == 'plasmid' else ''),
        xaxis_title="",
        xaxis=dict(
            tickmode='array',
            tickvals=[x_center],
            ticktext=[""]
        ),
        margin=dict(t=30, b=10, l=30, r=10),
        width=400,
        height=400,
        showlegend=False
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

# %%
def plot_gene(mge, gene, color, df_sub):
    df_sub[f'{gene}_count'] = df_sub[gene].apply(
        lambda x: 0 if pd.isna(x) or x.strip() == '' else len(x.split(','))
    )

    df_gene = df_sub.groupby(f'{gene}_count').size().reset_index(name='Count')

    x_max = max(1, df_gene[f'{gene}_count'].max())
    y_max = df_gene['Count'].max()
    y_dtick = tick_step(y_max)
    bar_width = 1

    fig = go.Figure()

    fig.add_trace(go.Bar(
        x=df_gene[f'{gene}_count'],
        y=df_gene['Count'],
        text=df_gene['Count'],
        textposition='inside',
        marker_color=color,
        width=[bar_width] * len(df_gene)
    ))
    
    fig.update_layout(
        title={
            'text': f'{gene} genes per {mge}',
            'y':1,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
            },
        xaxis=dict(
            type='linear',
            tickmode='linear',
            range=[-.5, x_max + .5],
            dtick=1,
            title=f'Number of {gene} genes'
        ),
        yaxis=dict(
            tickmode='linear',
            dtick=y_dtick,
            title=f'Number of {mge}s'
        ),
        margin=dict(t=30, b=10, l=30, r=10),
        width=400,
        height=400
    )
    
    return fig.to_html(full_html=False, include_plotlyjs='cdn')

# %%
def tick_step(max_val, max_ticks=10):
    """
    Compute a 'nice' integer tick step for y-axis,
    limiting number of ticks to max_ticks.
    """
    if max_val <= max_ticks:
        return 1
    
    rough_step = math.ceil(max_val / max_ticks)
    
    nice_steps = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
    for step in nice_steps:
        if rough_step <= step:
            return step
    return rough_step

# %%
def main():
    parser = argparse.ArgumentParser(
        description="Generate an interactive HTML report of MGEs from summary tables."
    )
    parser.add_argument(
        "-i", "--input-dir",
        required=True,
        help="Input directory with subfolders per sample, each containing MGE .tsv summaries and PNGs."
    )
    parser.add_argument(
        "-o", "--output-file",
        required=True,
        help="Output HTML report path."
    )
    args = parser.parse_args()

    input_dir = args.input_dir
    output_file = args.output_file

    # Scan dir, look for general reports
    report_files = [report for report in glob.glob(f'{input_dir}/*/*_summary.tsv')]
    df_complete = load_concat(report_files)

    # General dataframe (removing IS)
    df_general = df_complete[df_complete['MGE'] != 'IS']

    # Scan dir, look for specific MGE reports
    ## Plasmids
    plasmid_files = [report for report in glob.glob(f'{input_dir}/*/plasmids/*.tsv') if "plasmid_summary.tsv" not in report]
    df_pl = load_concat(plasmid_files)
    df_pl.rename(columns={
        'sample': 'Sample', 
        'contig': 'Contig',
        'mobsuite_id': 'Mobsuite ID',
        'ptu': 'PTU',
        'size': 'Length',
        'rep_mobsuite': 'Rep (MobS)',
        'rep_copla': 'Rep (Copla)',
        'mob_mobsuite': 'MOB (MobS)',
        'mob_copla': 'MOB (Copla)',
        'mpf_mobsuite': 'MPF (MobS)',
        'mpf_copla': 'MPF (Copla)',
        'orit_mobsuite': 'OriT',
        'predicted_mobility_mobsuite': 'Mobility (MobS)',
        'predicted_mobility_copla': 'Mobility (Copla)'
        }, inplace=True
    )

    ## Integrons
    integron_files = [report for report in glob.glob(f'{input_dir}/*/integrons/integron_summary.tsv')]
    df_int = load_concat(integron_files)

    ## Prophages
    phage_files = [report for report in glob.glob(f'{input_dir}/*/prophages/*.tsv') if "prophage_summary.tsv" not in report]
    df_ph = load_concat(phage_files)

    ## IS
    df_is = df_complete[df_complete['MGE'] == 'IS']
    df_is.drop(columns=["MGE", 'AMR', "VF", "DF"], inplace=True)
    df_is['Start'] = df_is['Start'].astype(int)
    df_is['End'] = df_is['End'].astype(int)

    # %%
    d_figs = {}

    d_figs['general'] = plot_general(df_general)

    df_filtered = df_general[
        (df_general['MGE'] != 'IS') &
        (~df_general['Name'].str.startswith('in0_', na=False))
    ].copy()

    for mge in df_filtered['MGE'].unique():
        df_sub = df_filtered[df_filtered['MGE'] == mge].copy()

        # Compute log length if plasmid, else use original
        df_sub['Length_log'] = (
            np.log10(df_sub['Length']) if mge == 'plasmid' else df_sub['Length']
        )

        d_figs[mge] = (plot_size(mge, df_sub), plot_gene(mge, 'AMR', '#FCAB10', df_sub), plot_gene(mge, 'VF', '#F8333C', df_sub), plot_gene(mge, 'DF', '#44AF69', df_sub))

    # %%
    render_final_report(
        df_general=df_filtered,
        df_pl=df_pl,
        df_int=df_int,
        df_ph=df_ph,
        df_is=df_is,
        d_figs=d_figs,
        input_dir=input_dir,
        output_file=output_file
    )

if __name__ == "__main__":
    main()
