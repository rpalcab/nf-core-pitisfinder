#!/usr/bin/env python3

import argparse
import pandas as pd
import glob
import plotly.express as px
import plotly.graph_objects as go
from jinja2 import Template
from datetime import datetime
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

{% macro render_section(title, mge, table_html, table_id, figs=[]) %}
<details class="section" open id="section-{{ mge }}">
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
<details class="section" open id="section-{{ mge }}">
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
        html {
            scroll-behavior: smooth;
        }
        #report-header {
            display: flex;
            align-items: center;
            background-color: #ffffff;
            border-bottom: 1px solid #dee2e6;
            padding: 1.2rem 2rem;
            margin-bottom: 2rem;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            font-family: 'Inter', sans-serif;
        }

        .header-left img {
            height: 120px;
            margin-right: 1.5rem;
        }

        .header-right h1 {
            margin: 0;
            font-size: 1.8rem;
            color: #333;
        }

        .header-right .subtitle {
            margin: 0.4rem 0;
            font-size: 1rem;
            color: #555;
            font-weight: 400;
        }

        .header-right .timestamp {
            margin-top: 0.5rem;
            font-size: 0.85rem;
            font-family: monospace;
            color: #888;
        }

        body {
            font-family: 'Inter', sans-serif;
            background-color: #ffffff;
            color: #212529;
            margin: 2rem;
            line-height: 1.6;
        }
        details.section {
            border: 1px solid #ddd;
            border-radius: 8px;
            padding: 1rem;
            margin-bottom: 2rem;
            transition: box-shadow 0.2s ease;
        }

        details.section:hover {
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
        }

        details > summary {
            cursor: pointer;
            font-size: 1.2rem;
            font-weight: 500;
            color: #333;
            list-style: none;
        }

        details[open] > summary {
            margin-bottom: 0.5rem;
        }
        .section {
            margin-bottom: 3rem;
        }
        details.section summary::before {
            content: "▶";
            display: inline-block;
            margin-right: 8px;
            transition: transform 0.2s ease;
        }

        details.section[open] summary::before {
            content: "▼";
        }
        table {
            margin-top: 1rem;
        }
        h2 {
            font-size: 1.4rem;
            margin: 0;
            color: #343a40;
        }
        .dataTables_wrapper { width: 100%; overflow-x: auto; }
        table.dataTable {
            table-layout: auto;
            background-color: #ffffff;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.05);
            overflow: hidden;
        }
        table.dataTable tbody tr:hover {
            background-color: #ffffff;
        }
        table.dataTable th,
        table.dataTable td {
            background-color: #ffffff;
            font-weight: 600;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            max-width: 240px;
            text-align: center !important;
        }
        table.dataTable td {
            font-weight: 400;
        }
        thead tr.filter-row th {
            padding: 4px;
            text-align: center;
        }
        .filter-box,
        .min-max-wrapper input,
        select.filter-box {
            width: 100%;
            box-sizing: border-box;
            padding: 0.4rem 0.6rem;
            border: 1px solid #ced4da;
            border-radius: 6px;
            background-color: #ffffff;
            font-size: 0.95rem;
            text-align: center;
        }
        .min-max-wrapper { display: flex; gap: 4px; }
        .dt-left-toolbar {
            display: flex;
            align-items: center;
            flex-wrap: wrap;
            gap: 0.5rem;
        }

        .sample-filter-wrapper {
            display: flex;
            align-items: center;
            gap: 0.5rem;
            margin-right: 1rem; /* space between filter and button */
        }
        input[type='number']::-webkit-outer-spin-button, input[type='number']::-webkit-inner-spin-button { -webkit-appearance: none; margin: 0; }
        input[type='number'] { -moz-appearance: textfield; }
        .plot-row {
            display: flex;
            flex-wrap: wrap;
            gap: 2rem;
            justify-content: center;
        }
        .plot-row iframe, .plotly-graph-div {
            max-width: 100%;
            height: auto;
        }
        details > summary {
            font-weight: 600;
            font-size: 1.2rem;
            padding: 0.5rem;
            cursor: pointer;
            list-style: none;
        }
        details > summary:hover { background-color: #ffffff; }
        details > summary h2 { display: inline; margin: 0; padding-left: 8px; font-size: 1.2em; }
        details > summary::-webkit-details-marker { display: none; }

        .fas.fa-image {
            color: #6c757d;
            transition: transform 0.2s ease, color 0.2s;
        }
        .fas.fa-image:hover {
            transform: scale(1.2);
            color: #343a40;
        }
        #floating-summary {
            position: fixed;
            top: 40px;
            right: 40px;
            background: #ffffffee;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            padding: 0.5rem 1rem;
            box-shadow: 0 2px 10px rgba(0,0,0,0.05);
            z-index: 999;
            font-size: 0.9rem;
            width: 150px;
        }

        #floating-summary summary {
            cursor: pointer;
            list-style: none;
            font-size: 1rem;
            color: #333;
        }

        #floating-summary ul {
            list-style: none;
            padding-left: 0;
            margin-top: 0.5rem;
        }

        #floating-summary ul li {
            margin: 6px 0;
        }

        #floating-summary a {
            color: #007bff;
            text-decoration: none;
            transition: color 0.2s;
        }

        @media (max-width: 900px) {
            #floating-summary {
                position: static;
                width: auto;
                margin-bottom: 2rem;
                box-shadow: none;
                border: none;
                background: transparent;
                padding: 0;
            }
        }
        .custom-select-sample {
            padding: 0.4rem 0.6rem;
            border: 1px solid #ced4da;
            border-radius: 6px;
            background-color: #ffffff;
            font-size: 0.95rem;
            color: #212529;
            font-family: 'Inter', sans-serif;
            appearance: none;
            cursor: pointer;
        }
        .custom-select-sample:focus {
            outline: none;
            border-color: #adb5bd;
            box-shadow: 0 0 0 0.1rem rgba(0,123,255,0.1);
        }
        div.dt-left-toolbar .dt-button {
            background-color: #ffffff !important;
            border: 1px solid #ced4da !important;
            border-radius: 6px !important;
            padding: 0.4rem 0.8rem !important;
            font-size: 0.95rem !important;
            font-weight: 500 !important;
            color: #212529 !important;
            font-family: 'Inter', sans-serif !important;
            cursor: pointer;
            transition: background-color 0.2s ease, box-shadow 0.2s ease;
        }

        div.dt-left-toolbar .dt-button:hover {
            background-color: #ffffff !important;
            color: #212529 !important;          /* enforce original text color */
            border-color: #ced4da !important;   /* enforce original border color */
            box-shadow: none !important;        /* remove any shadow on hover */
            outline: none !important;            /* remove focus outline */
        }
        .dt-left-toolbar label {
            font-size: 0.95rem;
            font-weight: 500;
            color: #343a40;
            margin-right: 0.5rem;
        }
        .custom-dt-button {
            background-color: #ffffff;
            margin-right: 0.5rem;
        }
        div.dt-button-collection {
            font-family: 'Inter', sans-serif;
            font-size: 0.95rem;
            border-radius: 6px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            border: 1px solid #dee2e6;
        }
        div.dt-button-collection button.dt-button {
            background-color: #ffffff;
            color: #212529;
            border: none;
            padding: 0.4rem 0.8rem;
            text-align: left;
        }

        div.dt-button-collection button.dt-button:hover {
            background-color: #ffffff !important;
        }

    </style>
</head>
<body>

<header id="report-header">
  <div class="header-left">
    <img src="https://raw.githubusercontent.com/rpalcab/nf-core-pitisfinder/main/docs/images/nf-core-pitisfinder_logo_light.png" alt="pitisfinder logo">
  </div>
  <div class="header-right">
    <h1>PITIsFinder (v.1.0)</h1>
    <p class="subtitle">Detection, characterization and classification of Mobile Genetic Elements (MGEs) from bacterial whole-genome assemblies.</p>
    <p class="timestamp">Report generated: {{ generation_date }}</p>
  </div>
</header>

<details id="floating-summary" open >
  <summary><strong>☰ Sections</strong></summary>
  <ul>
    <li><a href="#section-general">1. General Report</a></li>
    <li><a href="#section-plasmids">2. Plasmids</a></li>
    <li><a href="#section-integrons">3. Integrons</a></li>
    <li><a href="#section-prophages">4. Prophages</a></li>
    <li><a href="#section-IS">5. IS</a></li>
    <li><a href="#section-samples">6. Sample overview</a></li>
  </ul>
</details>

{{ render_section("1. General Report", 'general', render_table(df_general, zipped_columns_general, filter_options_general, 'table_general'), 'table_general', [d_figs.get('general')]) }}

{% if df_pl is not none %}
    {{ render_section("2. Plasmids", 'plasmids', render_table(df_pl, zipped_columns_pl, filter_options_pl, 'table_pl'), 'table_pl', d_figs.get('plasmid', [])) }}
{% else %}
    {{ render_empty('2. Plasmids', 'plasmids') }}
{% endif %}

{% if df_int is not none %}
    {{ render_section("3. Integrons", 'integrons', render_table(df_int, zipped_columns_int, filter_options_int, 'table_int'), 'table_int', d_figs.get('integron', [])) }}
{% else %}
    {{ render_empty('3. Integrons', 'integrons') }}
{% endif %}

{% if df_ph is not none %}
    {{ render_section("4. Prophages", 'prophages', render_table(df_ph, zipped_columns_ph, filter_options_ph, 'table_ph'), 'table_ph', d_figs.get('prophage', [])) }}
{% else %}
    {{ render_empty('4. Prophages', 'prophages') }}
{% endif %}

{% if df_is is not none %}
    {{ render_section("5. IS", 'IS', render_table(df_is, zipped_columns_is, filter_options_is, 'table_is'), 'table_is', d_figs.get('IS', [])) }}
{% else %}
    {{ render_empty('5. IS', 'IS') }}
{% endif %}

<details class="section" open id="section-samples">
    <summary>
        <span class="chevron"></span>
        <h2>6. Sample overview</h2>
    </summary>
    <div id="sample-filter-container-sample" class="sample-filter-container" style="margin-top: 10px;"></div>
</details>

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
                    text: 'Shown columns',
                    className: 'custom-dt-button'
                }
            ],
            initComplete: function () {
                const api = this.api();

                // 1) Sample‑column dropdown
                const sampleColIndex = 0;
                const uniqueSamples = api.column(sampleColIndex).data().unique().sort();
                const select = $(
                '<select class="custom-select-sample">' +
                '<option value="">All Samples</option>' +
                '</select>'
                );
                const container = api.table().container();
                uniqueSamples.each(d => {
                    if (d) select.append(`<option value="${d}">${d}</option>`);
                });
                const $label = $('<label class="me-2 mb-0">Select sample:</label>');
                const $wrapper = $('<div class="sample-filter-wrapper d-flex align-items-center me-3"></div>');
                $wrapper.append($label).append(select);

                $(api.table().container()).find('.dt-left-toolbar').prepend($wrapper);
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
def render_final_report(df_general, df_pl, df_int, df_ph, df_is, d_figs, generation_date, output_file):
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
        generation_date=generation_date
    )

    with open(output_file, "w") as f:
        f.write(rendered)

# %%
def load_concat(reports):
    if len(reports) == 0:
        return pd.DataFrame(columns=['Sample', 'Contig', 'Name'])
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
def infer_copy_png(df, png_paths, outpath):
    for _, row in df.iterrows():
        png = next(f for f in png_paths if f.name.endswith(row['tmp_img']))
        shutil.copy2(png, outpath / row['tmp_img'])
        df.at[_, 'Preview'] = outpath / row['tmp_img']
    df.drop(columns=['tmp_img'], inplace=True)
    df = inject_preview_icons(df)
    return df
# %%
def add_png_path(df, png_paths, output_file):
    outpath = Path(output_file).parent / "images"
    outpath.mkdir(parents=True, exist_ok=True)

    # Copy images to output directory, add path to df
    df['tmp_img'] = df.apply(lambda row: f"{row['Name']}.png" if row['MGE'] != 'plasmid' else f'{row['Name'].split(':')[0]}_{row['Sample']}.png', axis=1)
    df = infer_copy_png(df, png_paths, outpath)
    return df

# %%
def main():
    parser = argparse.ArgumentParser(
        description="Generate an interactive HTML report of MGEs from summary tables."
    )
    # parser.add_argument(
    #     "-i", "--input-dir",
    #     required=True,
    #     help="Input directory with subfolders per sample, each containing MGE .tsv summaries and PNGs."
    # )
    parser.add_argument(
        "-o", "--output_file",
        required=True,
        help="Output HTML report path."
    )
    parser.add_argument(
        "-g", "--gral_tsv",
        required=True,
        help="Path to the general summary TSV files of the samples."
    )
    parser.add_argument(
        "-G", "--gral_png",
        required=True,
        help="Path to the general PNG files of the samples."
    )
    parser.add_argument(
        "-p", "--plasmid_tsv",
        required=False,
        default="",
        help="Path to the plasmid summary TSV files of the samples."
    )
    parser.add_argument(
        "-P", "--plasmid_png",
        required=False,
        default="",
        help="Path to the plasmid PNG files of the samples."
    )
    parser.add_argument(
        "--plasmid_reports",
        required=False,
        default="",
        help="Path to individual plasmid reports TSV files of the samples."
    )
    parser.add_argument(
        "-i", "--integron_tsv",
        required=False,
        default="",
        help="Path to the integron summary TSV files of the samples."
    )
    parser.add_argument(
        "-I", "--integron_png",
        required=False,
        default="",
        help="Path to the integron PNG files of the samples."
    )
    parser.add_argument(
        "-f", "--phage_tsv",
        required=False,
        default="",
        help="Path to the phage summary TSV files of the samples."
    )
    parser.add_argument(
        "-F", "--phage_png",
        required=False,
        default="",
        help="Path to the phage PNG files of the samples."
    )
    args = parser.parse_args()

    output_file = args.output_file

    ## Sample reports
    report_files = [f for f in args.gral_tsv.split(',')]
    df_complete = load_concat(report_files)

    # Contig plots (chr/plasmids or contig_N)
    contig_plots = [Path(f) for f in args.gral_png.split(',')]
    plasmid_plots = [Path(f) for f in args.plasmid_png.split(',') if f != ""]
    list_plots = contig_plots + plasmid_plots

    # General dataframe (removing IS)
    df_general = df_complete[df_complete['MGE'] != 'IS']

    # ## Plasmids
    plasmid_files = [f for f in args.plasmid_reports.split(',') if f != ""]
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
    integron_files = [Path(f) for f in args.integron_tsv.split(',') if f != ""]
    df_int = load_concat(integron_files)

    ## Prophages
    phage_files = [Path(f) for f in args.phage_tsv.split(',') if f != ""]
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
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    integron_plots = [Path(f) for f in args.integron_png.split(',') if f != ""]
    phage_plots = [Path(f) for f in args.phage_png.split(',') if f != ""]
    mge_plots = plasmid_plots + integron_plots + phage_plots
    df_filtered = add_png_path(df_filtered, mge_plots, output_file)

    render_final_report(
        df_general=df_filtered,
        df_pl=df_pl,
        df_int=df_int,
        df_ph=df_ph,
        df_is=df_is,
        d_figs=d_figs,
        generation_date=now,
    #     input_dir=input_dir,
        output_file=output_file
    )

if __name__ == "__main__":
    main()
