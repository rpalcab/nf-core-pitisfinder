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
        'attI': "#F58230",
        'attC': "#F032E6",

        # Prophages (someday...)
        #TODO: When functional annotation is implemented, uncomment
        'Integrase': "#3CB44B",
        'Phage protein': "#F58230"
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
    gv = GenomeViz(fig_track_height=0.7, feature_track_ratio=0.5, track_align_type="center")
    gv.set_scale_bar(ymargin=0.5)

    track = gv.add_feature_track(gbk.records[0].id, gbk.get_seqid2size(), align_label=False)
    track.add_subtrack(name="GCcontent", ylim=(0, 100))
    track.add_subtrack(name="GCskew", ylim=(-1, 1))
    track.add_sublabel()

    # Track which legend elements are used
    legend_elements = set()

    # MGE composition
    features = gbk.extract_features(feature_type=['CDS', 'regulatory'])
    for feature in features:
        if feature.qualifiers.get('mge_element', [""])[0] == "yes":
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
                track.add_features(feature, hatch='///', fc="#0082C8", label_type="gene", ls="none")
                legend_elements.add('AMR')
            if 'VF' in feature.qualifiers.get('tag', [''])[0]:
                track.add_features(feature, fc="#CF0B0B", label_type="gene", ls="none")
                legend_elements.add('VF')
            elif 'hypothetical' in feature.qualifiers.get('product', [''])[0].lower():
                track.add_features(feature, label_type="gene", fc='grey', ls="none")
                legend_elements.add('hypothetical')
            else:
                track.add_features(feature, label_type="gene", fc="#0082C8", ls="none")
                legend_elements.add('CDS')

    # GCskew and content
    fig = gv.plotfig()
    win_size = 100
    step_size = 50
    gc_content_subtrack = track.get_subtrack("GCcontent")
    gc_skew_subtrack = track.get_subtrack("GCskew")
    seq = gbk.get_seqid2seq()[gbk.records[0].id]
    x, gc_content = gbk.calc_gc_content(window_size=win_size, step_size=step_size, seq=seq)
    x = track.segments[0].transform_coord(x)
    gc_content_subtrack.ax.fill_between(x, gc_content, color="grey")
    x, gc_skew = gbk.calc_gc_skew(window_size=win_size, step_size=step_size, seq=seq)
    x = track.segments[0].transform_coord(x)
    gc_skew_subtrack.ax.fill_between(x, gc_skew, color="pink")

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
        'AMR': ("#0082C8", '///', 'AMR gene'),
        'VF': ("#CF0B0B", None, 'VF gene'),
        'hypothetical': ("grey", None, 'Hypothetical gene'),
    }

    # Create handles only for present elements in consistent order
    ordered_keys = [
        'intI', 'Pc', 'Pint', 'attI', 'attC',
        #TODO: When functional annotation is implemented, modify to match new markers
        'Integrase', 'Phage protein',
        'AMR', 'VF', 'CDS', 'hypothetical'
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
