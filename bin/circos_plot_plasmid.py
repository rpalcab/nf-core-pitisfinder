#!/usr/bin/env python3

"""
Script to plot circular genomes based on GBK annotations.
Part of PITISfinder. Processes /tag qualifiers in features (if any)
"""

from pycirclize import Circos
from pycirclize.parser import Genbank
from Bio.SeqFeature import SeqFeature, FeatureLocation
import numpy as np
import pandas as pd
import math
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import argparse
from pathlib import Path

def get_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        prog = 'circos_plot_plasmid.py',
        description = 'circos_plot_plasmid.py is part of PITISfinder.'
    )

    parser.add_argument(
        '-i', '--input', required=True,
        type=Path,
        help="Required. GBK file"
    )

    parser.add_argument(
        '-m', '--mge', required=True,
        type=Path,
        help='Required. MGE coordinates table'
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

    df_mge = pd.read_csv(args.mge, sep="\t", header=0)

    gbk = Genbank(args.input)
    seqid2size = gbk.get_seqid2size()
    genome_size = sum(seqid2size.values())
    figsize, dpi = get_plot_params(genome_size)
    space = 0 if len(seqid2size) == 1 else 2
    circos = Circos(sectors=seqid2size, space=space)
    circos.text(f"{gbk.name}", size=12, r=125)

    seqid2features = gbk.get_seqid2features(feature_type=None)

    feature_presence = {
        "Forward CDS": False,
        "Reverse CDS": False,
        "Integron": False,
        "Prophage": False,
        "IS": False,
        "AMR": False,
        "VF": False,
        "DF": False,
        "MPF": False,
        "oriT": False,
        "MOB": False,
        "Replicon": False
    }

    track_info = [
        ("cds_track", (92, 98), "CDS"),
        ("int_track", (89, 92), "Integron"),
        ("ph_track", (86, 89), "Prophage"),
        ("is_track", (83, 86), "IS"),
        ("rvd_track", (80, 83), "AMR/VF/DF")
    ]

    mge_colors = {
        "MPF": "#9467BD",
        "oriT": "#FFBB78",
        "MOB": "#0F3B5A",
        "Replicon": "#912A2A"
    }

    for sector in circos.sectors:
        #Outer track
        outer_track = sector.add_track((98, 100))
        outer_track.axis(fc="lightgrey")

        # MGE track
        marker_track = sector.add_track((77, 80), r_pad_ratio=0.1)
        marker_track.axis(fc="#eaeaea", ec="lightgrey", lw=0.3)

        tracks = {'outer_track': outer_track}


        for i, (track_var, radius_range, label) in enumerate(track_info):
            bg_color = "#f8f8f8" if (i % 2 == 0) else "#eaeaea"
            track = sector.add_track(radius_range, r_pad_ratio=0.1)
            track.axis(fc=bg_color, ec="lightgrey", lw=0.3)
            tracks[track_var] = track

        features = seqid2features[sector.name]

        # Plot xticks & intervals on inner position
        contig_size = sector.size
        if contig_size > 1000000:
            div=1000000
            units="Mb"
        else:
            div=1000
            units="Kb"

        num_ticks = 13
        x = [int(i * contig_size / (num_ticks - 1)) for i in range(num_ticks)][:-1]
        labels = [f"{pos / div:.1f} {units}" for pos in x]

        marker_track.xticks(
            x=x,
            labels=labels,
            outer=False,
            label_orientation="vertical",
            line_kws=dict(ec="grey")
        )

        # Plot 'gene' qualifier label if exists
        labels, label_pos_list = [], []
        for feature in features:
            gene_name = feature.qualifiers.get("gene", [None])[0]
            if gene_name is not None and any(tag in feature.qualifiers.get('tag', []) for tag in ['AMR', 'VF', 'DF']):
                label_pos = (int(feature.location.start) + int(feature.location.end)) // 2
                tracks["cds_track"].annotate(label_pos, gene_name, label_size=7)
                labels.append(gene_name)
                label_pos_list.append(label_pos)

        # Plot CDS (fwd, rev) rRNA, tRNA and MGEs

        for feature in features:
            mge_subset = df_mge[
                (df_mge["MGE"] != "plasmid") &
                (df_mge["Contig"] == sector.name)
            ]

            for _, row in mge_subset.iterrows():
                start = int(row["Start"])
                end = int(row["End"])
                label = str(row["Name"])
                mge_type = row["MGE"].lower()

                # Create a mock SeqFeature
                fake_feature = SeqFeature(
                    FeatureLocation(start, end),
                    type="MGE",
                    qualifiers={"label": [label], "type": [mge_type]}
                )

                if mge_type == "integron":
                    track = tracks["int_track"]
                    color = "#FF7F0E"
                    feature_presence["Integron"] = True
                elif mge_type == "prophage":
                    track = tracks["ph_track"]
                    color = "#17BECF"
                    feature_presence["Prophage"] = True
                elif mge_type == "is":
                    track = tracks["is_track"]
                    color = "#E27FE4"
                    feature_presence["IS"] = True
                else:
                    continue

                # Plot the fake SeqFeature
                track.genomic_features([fake_feature], color=color, lw=0.1)

                # Annotate the name on CDS track
                label_pos = (start + end) // 2
                tracks["cds_track"].annotate(label_pos, label, label_size=7)

            # Gral features
            if feature.type == "CDS" and feature.location.strand == 1:
                tracks["cds_track"].genomic_features(feature, plotstyle="arrow", fc="#0082C8")
                feature_presence["Forward CDS"] = True
            elif feature.type == "CDS" and feature.location.strand == -1:
                tracks["cds_track"].genomic_features(feature, plotstyle="arrow", fc="#E6194B")
                feature_presence["Reverse CDS"] = True

            # AMR, VF, DF
            if feature.type == "CDS" and 'AMR' in feature.qualifiers.get('tag', []):
                tracks["rvd_track"].genomic_features(feature, color="#2CA02C", lw=0.1)
                feature_presence["AMR"] = True
            if feature.type == "CDS" and 'VF' in feature.qualifiers.get('tag', []):
                tracks["rvd_track"].genomic_features(feature, color="#FFD700", lw=0.1)
                feature_presence["VF"] = True
            if feature.type == "CDS" and 'DF' in feature.qualifiers.get('tag', []):
                tracks["rvd_track"].genomic_features(feature, color="#7F7F7F", lw=0.1)
                feature_presence["DF"] = True

            if feature.type == "CDS" and 'yes' in feature.qualifiers.get('mge_element', []):
                tag = feature.qualifiers['tag'][0]
                marker_track.genomic_features(feature, color=mge_colors[tag], lw=0.1)
                feature_presence[tag] = True
                start, end = int(feature.location.start), int(feature.location.end)
                label_pos = (start + end) // 2
                gene_name = feature.qualifiers.get("gene", [None])[0]
                tracks["cds_track"].annotate(label_pos, gene_name, label_size=7)

    fig = circos.plotfig(figsize=figsize)
    # Add legend
    handles = []
    legend_map = {
        "Forward CDS": Patch(color="#0082C8", label="Forward CDS"),
        "Reverse CDS": Patch(color="#E6194B", label="Reverse CDS"),
        "Integron": Patch(color="#FF7F0E", label="Integron"),
        "Prophage": Patch(color="#17BECF", label="Prophage"),
        "IS": Patch(color="#E27FE4", label="IS"),
        "AMR": Patch(color="#2CA02C", label="AMR"),
        "VF": Patch(color="#FFD700", label="Virulence Factor"),
        "DF": Patch(color="#7F7F7F", label="Defense Factor"),
        "MPF": Patch(color="#9467BD", label="MPF"),
        "oriT": Patch(color="#FFBB78", label="oriT"),
        "MOB": Patch(color="#0F3B5A", label="MOB"),
        "Replicon": Patch(color="#912A2A", label="Replicon")
    }

    for key, entry in legend_map.items():
        if feature_presence[key]:
            handles.append(entry)

    circos.ax.legend(handles=handles, bbox_to_anchor=(0.5, 0.5), loc="center", fontsize=12)

    fig.savefig(args.output, dpi=dpi)

if __name__ == "__main__":
    main()
