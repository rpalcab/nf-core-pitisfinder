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
        prog = 'circos_plot.py',
        description = 'linear_plot.py is part of PITISfinder.'
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

    #Data setup
    gv = GenomeViz(fig_track_height=0.7, feature_track_ratio=0.5, track_align_type="center")
    gv.set_scale_bar(ymargin=0.5)

    track = gv.add_feature_track(gbk.name, gbk.get_seqid2size(), align_label=False)
    track.add_subtrack(name="GCcontent", ylim=(0, 100))
    track.add_subtrack(name="GCskew", ylim=(-1, 1))
    track.add_sublabel()

    # EGM composition
    features = gbk.extract_features(feature_type=['CDS', 'regulatory'])
    for feature in features:
        if feature.type == "CDS" and 'integrase' in feature.qualifiers.get('product', [''])[0].lower():
            track.add_features(feature, label_type="gene", fc='blue', lw=1.0)
        elif feature.type == "CDS" and 'AMR' in feature.qualifiers.get('tag', [''])[0]:
            track.add_features(feature, hatch='//', fc='orange', label_type="gene", lw=1.0)
        elif feature.type == "CDS" and 'hypothetical' in feature.qualifiers.get('product', [''])[0].lower():
            track.add_features(feature, label_type="gene", fc='grey', lw=1.0)
        elif feature.type == "CDS":
            track.add_features(feature, label_type="gene", fc='orange', lw=1.0)
        elif feature.type == "regulatory":
            if 'atti' in feature.qualifiers.get('note', [''])[0].lower():
                feature.qualifiers['note'] = ['attI']
            track.add_features(feature, label_type="note", fc='green', lw=1.0)

    # GCskew and content
    fig = gv.plotfig()
    win_size = 100
    step_size = 50
    gc_content_subtrack = track.get_subtrack("GCcontent")
    gc_skew_subtrack = track.get_subtrack("GCskew")
    seq = gbk.get_seqid2seq()[gbk.name]
    x, gc_content = gbk.calc_gc_content(window_size=win_size, step_size=step_size, seq=seq)
    x = track.segments[0].transform_coord(x)
    gc_content_subtrack.ax.fill_between(x, gc_content, color="grey")
    x, gc_skew = gbk.calc_gc_skew(window_size=win_size, step_size=step_size, seq=seq)
    x = track.segments[0].transform_coord(x)
    gc_skew_subtrack.ax.fill_between(x, gc_skew, color="pink")

    # Legend
    _ = fig.legend(
        handles=[
            Patch(facecolor="blue", label="Integrase"),
            Patch(facecolor="orange", label="CDS"),
            Patch(facecolor="orange", hatch="///", label="AMR gene"),
            Patch(facecolor="grey", label="Hypothetical gene"),
            Patch(facecolor="green", label="att site")
        ],
        fontsize=12,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        handlelength=1.0,
    )

    fig.savefig(args.output)

if __name__ == "__main__":
    main()
