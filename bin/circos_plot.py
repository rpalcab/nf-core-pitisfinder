#!/usr/bin/env python3

"""
Script to plot circular genomes based on GBK annotations.
Part of PITISfinder. Processes /tag qualifiers in features (if any)
"""

from pycirclize import Circos
from pycirclize.parser import Genbank
import numpy as np
import math
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
        '-m', '--mge_elements',
        action="store_true",
        help="Plot MGE elements (integrase, MPF, MOB...)"
    )

    parser.add_argument(
        '-o', '--output', required=True,
        type=Path,
        help='Required. Final output png file'
    )

    return parser.parse_args()

# %%
def get_plot_params(genome_size):
    # Cap sizes to avoid too-small or too-huge figures
    base_size = 12
    scale_factor = np.log10(genome_size / 100000)  # base 6 for 100 kb
    scale_factor = max(0, min(scale_factor, 5))  # Clamp between 0 and 5
    figsize = (base_size + scale_factor, base_size + scale_factor)
    dpi = int(300 + 50 * scale_factor)
    return tuple(figsize), dpi

# %%
def main():
    args = get_args()

    gbk = Genbank(args.input)
    seqid2size = gbk.get_seqid2size()
    genome_size = sum(seqid2size.values())
    figsize, dpi = get_plot_params(genome_size)
    space = 0 if len(seqid2size) == 1 else 2
    circos = Circos(sectors=seqid2size, space=space)
    circos.text(f"{gbk.name}", size=12, r=150)

    seqid2features = gbk.get_seqid2features(feature_type=None)
    seqid2seq = gbk.get_seqid2seq()

    feature_presence = {
        "Forward CDS": False,
        "Reverse CDS": False,
        "rRNA": False,
        "tRNA": False,
        "Plasmid": False,
        "Integron": False,
        "Prophage": False,
        "IS": False,
        "AMR": False,
        "VF": False,
        "DF": False,
        "MPF": False,
        "oriT": False,
        "MOB": False,
        "Replicon": False,
        "GC Content": True,
        "GC Skew": True,
    }

    track_info = [
        ("cds_track", (92, 98), "CDS"),
        ("rna_track", (89, 92), "rRNA/tRNA"),
        ("pl_track", (86, 89), "Plasmid"),
        ("int_track", (83, 86), "Integron"),
        ("ph_track", (80, 83), "Prophage"),
        ("is_track", (77, 80), "IS"),
        ("rvd_track", (74, 77), "AMR/VF/DF")
    ]

    mge_colors = {
        "MPF": "#911EB4",
        "oriT": "#F58230",
        "MOB": "#46F0F0",
        "Replicon": "#F032E6"
    }

    for sector in circos.sectors:
        #Outer track
        outer_track = sector.add_track((98, 100))
        outer_track.axis(fc="lightgrey")
        # GC tracks
        gc_content_track = sector.add_track((65, 70), r_pad_ratio=0.1)
        gc_skew_track = sector.add_track((55, 65), r_pad_ratio=0.1)

        # MGE track (if wanted)
        if args.mge_elements is True:
            marker_track = sector.add_track((71, 74), r_pad_ratio=0.1)
            marker_track.axis(fc="#eaeaea", ec="lightgrey", lw=0.3)

        tracks = {'outer_track': outer_track,
                  'gc_content_track': gc_content_track,
                  'gc_skew_track': gc_skew_track                  }


        for i, (track_var, radius_range, label) in enumerate(track_info):
            bg_color = "#f8f8f8" if (i % 2 == 0) else "#eaeaea"
            track = sector.add_track(radius_range, r_pad_ratio=0.1)
            track.axis(fc=bg_color, ec="lightgrey", lw=0.3)
            tracks[track_var] = track

        features = seqid2features[sector.name]

        # Plot xticks & intervals on inner position
        genome_size = sum(gbk.get_seqid2size().values())
        if genome_size > 1000000:
            interval = math.ceil((genome_size * .97 // 8) / 100000) * 100000
            div=1000000
            units="Mb"
        else:
            interval = math.ceil((genome_size // 8) / 1000) * 1000
            div=1000
            units="Kb"
        tracks["cds_track"].xticks_by_interval(
            interval=interval,
            outer=True,
            label_formatter=lambda v: f"{v/ div:.1f} {units}",
            label_orientation="vertical",
            line_kws=dict(ec="grey"),
        )

        # Plot 'gene' qualifier label if exists
        labels, label_pos_list = [], []
        for feature in features:
            start, end = int(feature.location.start), int(feature.location.end)
            label_pos = (start + end) // 2
            gene_name = feature.qualifiers.get("gene", [None])[0]
            if gene_name is not None and any(tag in feature.qualifiers.get('tag', []) for tag in ['AMR', 'VF', 'DF']):
                tracks["cds_track"].annotate(label_pos, gene_name, label_size=7, text_kws = {"weight":"bold"})
                labels.append(gene_name)
                label_pos_list.append(label_pos)

        # Plot CDS (fwd, rev) rRNA, tRNA and MGEs

        for feature in features:
            # Gral features
            if feature.type == "CDS" and feature.location.strand == 1:
                tracks["cds_track"].genomic_features(feature, plotstyle="arrow", fc="#0082C8")
                feature_presence["Forward CDS"] = True
            elif feature.type == "CDS" and feature.location.strand == -1:
                tracks["cds_track"].genomic_features(feature, plotstyle="arrow", fc="#E6194B")
                feature_presence["Reverse CDS"] = True
            elif feature.type == "rRNA":
                tracks["rna_track"].genomic_features(feature, fc="#3CB44B")
                feature_presence["rRNA"] = True
            elif feature.type == "tRNA":
                tracks["rna_track"].genomic_features(feature, color="#FFE119", lw=0.1)
                feature_presence["tRNA"] = True

            # MGEs
            elif feature.type == "MGE" and 'plasmid' in feature.qualifiers.get('type', []):
                tracks["pl_track"].genomic_features(feature, color="#911EB4", lw=0.1)
                feature_presence["Plasmid"] = True
            elif feature.type == "MGE" and 'integron' in feature.qualifiers.get('type', []):
                tracks["int_track"].genomic_features(feature, color="#F58230", lw=0.1)
                feature_presence["Integron"] = True
            elif feature.type == "MGE" and 'prophage' in feature.qualifiers.get('type', []):
                tracks["ph_track"].genomic_features(feature, color="#46F0F0", lw=0.1)
                feature_presence["Prophage"] = True
            elif feature.type == "MGE" and 'IS' in feature.qualifiers.get('type', []):
                tracks["is_track"].genomic_features(feature, color="#F032E6", lw=0.1)
                feature_presence["IS"] = True

            # AMR, VF, DF
            if feature.type == "CDS" and 'AMR' in feature.qualifiers.get('tag', []):
                tracks["rvd_track"].genomic_features(feature, color="#D2F53C", lw=0.1)
                feature_presence["AMR"] = True
            if feature.type == "CDS" and 'VF' in feature.qualifiers.get('tag', []):
                tracks["rvd_track"].genomic_features(feature, color="#AA6E28", lw=0.1)
                feature_presence["VF"] = True
            if feature.type == "CDS" and 'DF' in feature.qualifiers.get('tag', []):
                tracks["rvd_track"].genomic_features(feature, color="#808080", lw=0.1)
                feature_presence["DF"] = True

            # MGE biomarkers (if wanted)
            if args.mge_elements is True:
                if feature.type == "CDS" and 'yes' in feature.qualifiers.get('mge_element', []):
                    tag = feature.qualifiers['tag'][0]
                    marker_track.genomic_features(feature, color=mge_colors[tag], lw=0.1)
                    feature_presence[tag] = True
                    start, end = int(feature.location.start), int(feature.location.end)
                    label_pos = (start + end) / 2
                    gene_name = feature.qualifiers.get("gene", [None])[0]
                    tracks["cds_track"].annotate(label_pos, gene_name, label_size=7)

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

    fig = circos.plotfig(figsize=figsize)
    # Add legend
    handles = []

    legend_map = {
        "Forward CDS": Patch(color="#E6194B", label="Forward CDS"),
        "Reverse CDS": Patch(color="#0082C8", label="Reverse CDS"),
        "rRNA": Patch(color="#3CB44B", label="rRNA"),
        "tRNA": Patch(color="#FFE119", label="tRNA"),
        "Plasmid": Patch(color="#911EB4", label="Plasmid"),
        "Integron": Patch(color="#F58230", label="Integron"),
        "Prophage": Patch(color="#46F0F0", label="Prophage"),
        "IS": Patch(color="#F032E6", label="IS"),
        "AMR": Patch(color="#D2F53C", label="AMR"),
        "VF": Patch(color="#AA6E28", label="VF"),
        "DF": Patch(color="#808080", label="DF"),
        "MPF": Patch(color="#911EB4", label="MPF"),
        "oriT": Patch(color="#F58230", label="oriT"),
        "MOB": Patch(color="#46F0F0", label="MOB"),
        "Replicon": Patch(color="#F032E6", label="Replicon"),
        "GC Content": [
            Line2D([], [], color="black", label="Positive GC Content", marker="^", ms=6, ls="None"),
            Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=6, ls="None"),
        ],
        "GC Skew": [
            Line2D([], [], color="olive", label="Positive GC Skew", marker="^", ms=6, ls="None"),
            Line2D([], [], color="purple", label="Negative GC Skew", marker="v", ms=6, ls="None"),
        ],
    }

    for key, entry in legend_map.items():
        if feature_presence[key]:
            if isinstance(entry, list):
                handles.extend(entry)
            else:
                handles.append(entry)

    circos.ax.legend(handles=handles, bbox_to_anchor=(0.5, 0.5), loc="center", fontsize=8)

    fig.savefig(args.output, dpi=dpi)

if __name__ == "__main__":
    main()
