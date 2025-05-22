#!/usr/bin/env python3

"""
Script to plot circular genomes based on GBK annotations.
Part of PITISfinder. Processes /tag qualifiers in features (if any)
"""

from pycirclize import Circos
from pycirclize.parser import Genbank
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import argparse
from pathlib import Path

def get_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        prog = 'circos_plot.py',
        description = 'circos_plot.py is part of PITISfinder.'
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

# %%
def main():
    args = get_args()

    gbk = Genbank(args.input)
    seqid2size = gbk.get_seqid2size()
    space = 0 if len(seqid2size) == 1 else 2
    circos = Circos(sectors=seqid2size, space=space)
    circos.text(f"{gbk.name}", size=12, r=20)

    seqid2features = gbk.get_seqid2features(feature_type=None)
    seqid2seq = gbk.get_seqid2seq()
    for sector in circos.sectors:
        # Plot outer track with xticks
        major_ticks_interval = 500000
        minor_ticks_interval = 100000
        outer_track = sector.add_track((98, 100))
        outer_track.axis(fc="lightgrey")
        if sector.size >= major_ticks_interval:
            outer_track.xticks_by_interval(
                major_ticks_interval, label_formatter=lambda v: f"{v/ 10 ** 6:.1f} Mb"
            )
            outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)

        f_cds_track = sector.add_track((90, 97), r_pad_ratio=0.1)
        r_cds_track = sector.add_track((83, 90), r_pad_ratio=0.1)
        # rrna_track = sector.add_track((76, 83), r_pad_ratio=0.1)
        # trna_track = sector.add_track((76, 83), r_pad_ratio=0.1)
        amr_track = sector.add_track((76, 83), r_pad_ratio=0.1)
        pl_track = sector.add_track((73, 76), r_pad_ratio=0.1)
        int_track = sector.add_track((70, 73), r_pad_ratio=0.1)
        ph_track = sector.add_track((67, 70), r_pad_ratio=0.1)
        is_track = sector.add_track((64, 67), r_pad_ratio=0.1)

        # Plot Forward CDS, Reverse CDS, rRNA, tRNA
        features = seqid2features[sector.name]
        for feature in features:
            if feature.type == "CDS" and feature.location.strand == 1:
                f_cds_track.genomic_features(feature, fc="red")
            elif feature.type == "CDS" and feature.location.strand == -1:
                r_cds_track.genomic_features(feature, fc="blue")
            # elif feature.type == "rRNA":
            #     rrna_track.genomic_features(feature, fc="green")
            # if feature.type == "tRNA":
            #     trna_track.genomic_features(feature, color="magenta", lw=0.1)
            if feature.type == "CDS" and 'AMR' in feature.qualifiers.get('tag', []):
                amr_track.genomic_features(feature, color="cyan", lw=0.1)
            if feature.type == "MGE" and 'plasmid' in feature.qualifiers.get('type', []):
                pl_track.genomic_features(feature, color="green", lw=0.1)
            if feature.type == "MGE" and 'integron' in feature.qualifiers.get('type', []):
                int_track.genomic_features(feature, color="magenta", lw=0.1)
            if feature.type == "MGE" and 'prophage' in feature.qualifiers.get('type', []):
                ph_track.genomic_features(feature, color="grey", lw=0.1)
            if feature.type == "MGE" and 'IS' in feature.qualifiers.get('type', []):
                is_track.genomic_features(feature, color="purple", lw=0.1)

        # Plot 'gene' qualifier label if exists
        labels, label_pos_list = [], []
        for feature in features:
            start = int(feature.location.start)
            end = int(feature.location.end)
            label_pos = (start + end) / 2
            gene_name = feature.qualifiers.get("gene", [None])[0]
            if gene_name is not None and 'AMR' in feature.qualifiers.get('tag', []):
                labels.append(gene_name)
                label_pos_list.append(label_pos)
        f_cds_track.xticks(label_pos_list, labels, label_size=6, label_orientation="vertical")

        # # Plot xticks (interval = 10 Kb)
        # r_cds_track.xticks_by_interval(
        #     10000, outer=False, label_formatter=lambda v: f"{v/1000:.1f} Kb"
        # )

        # Plot GC content
        gc_content_track = sector.add_track((55, 60))
        seq = seqid2seq[sector.name]
        label_pos_list, gc_contents = gbk.calc_gc_content(seq=seq)
        gc_contents = gc_contents - gbk.calc_genome_gc_content(seq=gbk.full_genome_seq)
        positive_gc_contents = np.where(gc_contents > 0, gc_contents, 0)
        negative_gc_contents = np.where(gc_contents < 0, gc_contents, 0)
        abs_max_gc_content = np.max(np.abs(gc_contents))
        vmin, vmax = -abs_max_gc_content, abs_max_gc_content
        gc_content_track.fill_between(
            label_pos_list, positive_gc_contents, 0, vmin=vmin, vmax=vmax, color="black"
        )
        gc_content_track.fill_between(
            label_pos_list, negative_gc_contents, 0, vmin=vmin, vmax=vmax, color="grey"
        )

        # Plot GC skew
        gc_skew_track = sector.add_track((35, 50))

        label_pos_list, gc_skews = gbk.calc_gc_skew(seq=seq)
        positive_gc_skews = np.where(gc_skews > 0, gc_skews, 0)
        negative_gc_skews = np.where(gc_skews < 0, gc_skews, 0)
        abs_max_gc_skew = np.max(np.abs(gc_skews))
        vmin, vmax = -abs_max_gc_skew, abs_max_gc_skew
        gc_skew_track.fill_between(
            label_pos_list, positive_gc_skews, 0, vmin=vmin, vmax=vmax, color="olive"
        )
        gc_skew_track.fill_between(
            label_pos_list, negative_gc_skews, 0, vmin=vmin, vmax=vmax, color="purple"
        )

    fig = circos.plotfig()
    # Add legend
    handles = [
        Patch(color="red", label="Forward CDS"),
        Patch(color="blue", label="Reverse CDS"),
        # Patch(color="green", label="rRNA"),
        # Patch(color="magenta", label="tRNA"),
        Patch(color="cyan", label="AMR"),
        Patch(color="green", label="Plasmid"),
        Patch(color="magenta", label="Integron"),
        Patch(color="grey", label="Prophage"),
        Patch(color="purple", label="IS"),

        Line2D([], [], color="black", label="Positive GC Content", marker="^", ms=6, ls="None"),
        Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=6, ls="None"),
        Line2D([], [], color="olive", label="Positive GC Skew", marker="^", ms=6, ls="None"),
        Line2D([], [], color="purple", label="Negative GC Skew", marker="v", ms=6, ls="None"),
    ]
    _ = circos.ax.legend(handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=8)

    fig.savefig(args.output)

if __name__ == "__main__":
    main()
