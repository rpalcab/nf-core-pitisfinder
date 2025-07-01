#!/usr/bin/env python3
"""
Script to plot linear sequences based on GBK annotations.
Part of PITISfinder. Processes /tag qualifiers in features (if any)
"""

from pygenomeviz.parser import Genbank
from pathlib import Path
from pygenomeviz import GenomeViz
from matplotlib.patches import Patch
import argparse

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog='linear_plot.py',
        description='linear_plot.py is part of PITISfinder.'
    )
    parser.add_argument(
        '-i', '--input', required=True,
        type=Path,
        help="Required. GBK file"
    )
    parser.add_argument(
        '-o', '--output', required=True,
        type=Path,
        help='Required. Final output png file'
    )
    return parser.parse_args()

def main():
    args = get_args()
    gbk = Genbank(args.input)

    d_tag = {
        # Integrons
        'intI': "#3CB44B",
        'Pc': "#FFE119",
        'Pint': "#46F0F0",
        'attI': "#78608F",
        'attC': "#F032E6",

        # Prophages (someday...)
        #TODO: When functional annotation is implemented, uncomment
        'Integrase': "#3CB44B",
        'Phage protein': "#ECBD56"
        # 'Assembly': "#1F77B4",
        # 'Baseplate': "#AEC7E8",
        # 'Coat': "#FF7F0E",
        # 'Integration': "#9467BD",
        # 'Lysis': "#D62728",
        # 'Other': "#7F7F7F",
        # 'Other (structural)': "#BCBD22",
        # 'Portal': "#17BECF",
        # 'Replication': "#E377C2",
        # 'Tail': "#8C564B",
        # 'Terminase': "#2CA02C"
        # ICEs
    }

    # Data setup
    # Calculate dynamic plot width based on sequence length
    seq_len = gbk.get_seqid2size()[gbk.records[0].id]

    # Parameters you can tune
    base_width = 15       # Minimum width
    scale_factor = 0.0005 # Increase to reduce overlap (more space per bp)

    plot_width = max(base_width, seq_len * scale_factor)

    gv = GenomeViz(fig_width=plot_width, fig_track_height=0.7, feature_track_ratio=0.5, track_align_type="center")
    gv.set_scale_bar(ymargin=0.5)

    track = gv.add_feature_track(gbk.records[0].id, gbk.get_seqid2size(), align_label=False)
    track.add_sublabel()

    # Track which legend elements are used
    legend_elements = set()

    # MGE composition
    features = gbk.extract_features(feature_type=['CDS', 'regulatory'])
    for feature in features:
        if feature.qualifiers.get('mge_element', [""])[0] == "yes":
            product = feature.qualifiers.get('product', [''])[0]
            truncated_product = (product[:41] + '...') if len(product) > 50 else product
            feature.qualifiers['product'] = [truncated_product]
            tag = feature.qualifiers['tag'][0]
            if tag in d_tag and tag != "Phage protein":
                track.add_features(feature, label_type='tag', fc=d_tag[tag], ls="none")
                legend_elements.add(tag)
            #TODO: When functional annotation is implemented, remove this conditional block
            elif tag in d_tag and tag == "Phage protein":
                track.add_features(feature, label_type='product', fc=d_tag[tag], ls="none")
                legend_elements.add(tag)

        # General CDS
        elif feature.type == "CDS":
            if 'AMR' in feature.qualifiers.get('tag', [''])[0]:
                track.add_features(feature, fc="#F58230", label_type="gene", ls="none")
                legend_elements.add('AMR')
            elif 'VF' in feature.qualifiers.get('tag', [''])[0]:
                track.add_features(feature, fc="#CF0B0B", label_type="gene", ls="none")
                legend_elements.add('VF')
            elif 'DF' in feature.qualifiers.get('tag', [''])[0]:
                track.add_features(feature, fc="#6B0DA1", label_type="gene", ls="none")
                legend_elements.add('DF')
            elif 'hypothetical' in feature.qualifiers.get('product', [''])[0].lower():
                track.add_features(feature, label_type="gene", fc='grey', ls="none")
                legend_elements.add('hypothetical')
            else:
                track.add_features(feature, label_type="gene", fc="#0082C8", ls="none")
                legend_elements.add('CDS')

    fig = gv.plotfig()

    # Dynamic legend creation
    legend_map = {
        'intI': (d_tag['intI'], None, 'Integrase'),
        'Pc': (d_tag['Pc'], None, 'Promoter'),
        'Pint': (d_tag['Pint'], None, 'Integron promoter'),
        'attI': (d_tag['attI'], None, 'attI site'),
        'attC': (d_tag['attC'], None, 'attC site'),
        #TODO: When functional annotation is implemented, modify to match new markers
        'Integrase': (d_tag['Integrase'], None, 'Integrase'),
        'Phage protein': (d_tag['Phage protein'], None, 'Phage protein'),
        'CDS': ("#0082C8", None, 'CDS'),
        'AMR': ("#F58230", None, 'AMR gene'),
        'VF': ("#CF0B0B", None, 'VF gene'),
        'DF': ("#6B0DA1", None, 'DF gene'),
        'hypothetical': ("grey", None, 'Hypothetical gene'),
    }

    # Create handles only for present elements in consistent order
    ordered_keys = [
        'intI', 'Pc', 'Pint', 'attI', 'attC',
        #TODO: When functional annotation is implemented, modify to match new markers
        'Integrase', 'Phage protein',
        'AMR', 'VF', 'DF', 'CDS', 'hypothetical'
    ]
    legend_handles = []
    for key in ordered_keys:
        if key in legend_elements and key in legend_map:
            fc, hatch, label = legend_map[key]
            legend_handles.append(
                Patch(facecolor=fc, hatch=hatch, label=label)
            )

    # Add legend if any elements are present
    if legend_handles:
        fig.legend(
            handles=legend_handles,
            fontsize=12,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            handlelength=1.0,
        )

    fig.savefig(args.output, bbox_inches='tight')

if __name__ == "__main__":
    main()
