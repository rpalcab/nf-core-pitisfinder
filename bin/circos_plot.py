#!/usr/bin/env python3

"""
Script to plot circular genomes based on GBK annotations.
Part of PITISfinder. Processes /tag qualifiers in features (if any)
"""

from pycirclize import Circos
from pycirclize.parser import Genbank
import matplotlib as plt
import numpy as np
import pandas as pd
import math
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import argparse
from pathlib import Path

# %%
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
        "-r", "--mobsuite_report",
        default=None,
        help="MOBsuite report with plasmid and chromosome information"
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
def group_contigs_by_replicon(mobsuite_path):
    df = pd.read_csv(mobsuite_path, sep='\t')
    replicon_map = {}

    for _, row in df.iterrows():
        contig = row['contig_id']
        mol_type = row['molecule_type']
        cluster_id = row['primary_cluster_id']

        if mol_type == "chromosome":
            replicon_id = "chromosome"
        elif mol_type == "plasmid" and pd.notna(cluster_id) and cluster_id != "-":
            replicon_id = f"plasmid_{cluster_id}"
        else:
            continue  # Ignore unclustered plasmids

        replicon_map.setdefault(replicon_id, []).append(contig)

    return replicon_map

# %%
def handle_fragmented_assembly(gbk_parser, output_prefix):
    contig_sizes = gbk_parser.get_seqid2size()
    filtered_contigs = {k: v for k, v in contig_sizes.items() if v >= 1000}

    if len(filtered_contigs) > 15:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.text(0.5, 0.5, "Assembly too fragmented. Not plotted.",
                fontsize=16, ha='center', va='center')
        ax.axis('off')
        fig.savefig(f"{output_prefix}_too_fragmented.png", dpi=300)
        return

    for contig in filtered_contigs:
        contig_output = f"{output_prefix}.contig_{contig}.png"
        plot_single_contig(gbk_parser, contig, contig_output)

def plot_single_contig(gbk_parser, contig_id, output_path):
    single_parser = Genbank(gbk_parser.file_path, seq_ids=[contig_id])
    plot_circos(single_parser, output_path, mge_elements=False)

def plot_circos(seqid2size, seqid2features, title, output_path, mge_elements=False):
    genome_size = sum(seqid2size.values())
    figsize, dpi = get_plot_params(genome_size)
    space = 0 if len(seqid2size) == 1 else 2
    circos = Circos(sectors=seqid2size, space=space)
    circos.text(f"{title}", size=12, r=125)

    feature_presence = {
        "Forward CDS": False, "Reverse CDS": False, "tRNA": False,
        "Plasmid": False, "Integron": False, "Prophage": False, "IS": False,
        "AMR": False, "VF": False, "DF": False,
        "MPF": False, "oriT": False, "MOB": False, "Replicon": False
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
        "MPF": "#9467BD", "oriT": "#FFBB78",
        "MOB": "#0F3B5A", "Replicon": "#912A2A"
    }

    for sector in circos.sectors:
        outer_track = sector.add_track((98, 100))
        outer_track.axis(fc="lightgrey")
        marker_track = None

        if mge_elements:
            marker_track = sector.add_track((71, 74), r_pad_ratio=0.1)
            marker_track.axis(fc="#eaeaea", ec="lightgrey", lw=0.3)

        tracks = {'outer_track': outer_track}

        for i, (track_var, radius_range, label) in enumerate(track_info):
            bg_color = "#f8f8f8" if (i % 2 == 0) else "#eaeaea"
            track = sector.add_track(radius_range, r_pad_ratio=0.1)
            track.axis(fc=bg_color, ec="lightgrey", lw=0.3)
            tracks[track_var] = track

        features = seqid2features[sector.name]
        contig_size = sector.size
        div = 1_000_000 if contig_size > 1_000_000 else 1_000
        units = "Mb" if div == 1_000_000 else "Kb"

        num_ticks = 13
        x = [int(i * contig_size / (num_ticks - 1)) for i in range(num_ticks)][:-1]
        labels = [f"{pos / div:.1f} {units}" for pos in x]
        tracks["rvd_track"].xticks(x=x, labels=labels, outer=False,
                                   label_orientation="vertical", line_kws=dict(ec="grey"))

        for feature in features:
            start, end = int(feature.location.start), int(feature.location.end)
            label_pos = (start + end) // 2
            gene_name = feature.qualifiers.get("gene", [None])[0]

            # Annotate if AMR/VF/DF
            if gene_name and any(tag in feature.qualifiers.get('tag', []) for tag in ['AMR', 'VF', 'DF']):
                tracks["cds_track"].annotate(label_pos, gene_name, label_size=7, text_kws={"weight": "bold"})

        for feature in features:
            # CDS & RNA
            if feature.type == "CDS":
                color = "#0082C8" if feature.location.strand == 1 else "#E6194B"
                tracks["cds_track"].genomic_features(feature, plotstyle="arrow", fc=color)
                key = "Forward CDS" if feature.location.strand == 1 else "Reverse CDS"
                feature_presence[key] = True
            elif feature.type == "tRNA":
                tracks["rna_track"].genomic_features(feature, color="#10470B", lw=0.1)
                feature_presence["tRNA"] = True

            # MGE annotations
            if feature.type == "MGE":
                type_list = feature.qualifiers.get('type', [])
                if 'plasmid' in type_list:
                    tracks["pl_track"].genomic_features(feature, color="#911EB4", lw=0.1)
                    feature_presence["Plasmid"] = True
                elif 'integron' in type_list:
                    tracks["int_track"].genomic_features(feature, color="#FF7F0E", lw=0.1)
                    feature_presence["Integron"] = True
                elif 'prophage' in type_list:
                    tracks["ph_track"].genomic_features(feature, color="#17BECF", lw=0.1)
                    feature_presence["Prophage"] = True
                elif 'IS' in type_list:
                    tracks["is_track"].genomic_features(feature, color="#E27FE4", lw=0.1)
                    feature_presence["IS"] = True

            # Tags (AMR/VF/DF)
            for tag, color in zip(['AMR', 'VF', 'DF'], ['#2CA02C', '#FFD700', '#7F7F7F']):
                if feature.type == "CDS" and tag in feature.qualifiers.get('tag', []):
                    tracks["rvd_track"].genomic_features(feature, color=color, lw=0.1)
                    feature_presence[tag] = True

            # MGE markers
            if mge_elements and feature.type == "CDS" and 'yes' in feature.qualifiers.get('mge_element', []):
                tag = feature.qualifiers['tag'][0]
                marker_track.genomic_features(feature, color=mge_colors[tag], lw=0.1)
                feature_presence[tag] = True
                gene_name = feature.qualifiers.get("gene", [None])[0]
                if gene_name:
                    label_pos = (int(feature.location.start) + int(feature.location.end)) // 2
                    tracks["cds_track"].annotate(label_pos, gene_name, label_size=7)

    fig = circos.plotfig(figsize=figsize)
    # Build legend
    legend_map = {
        "Forward CDS": "#0082C8", "Reverse CDS": "#E6194B", "tRNA": "#10470B",
        "Plasmid": "#911EB4", "Integron": "#FF7F0E", "Prophage": "#17BECF", "IS": "#E27FE4",
        "AMR": "#2CA02C", "VF": "#FFD700", "DF": "#7F7F7F",
        "MPF": "#9467BD", "oriT": "#FFBB78", "MOB": "#0F3B5A", "Replicon": "#912A2A"
    }
    handles = [Patch(color=color, label=key) for key, color in legend_map.items() if feature_presence[key]]
    circos.ax.legend(handles=handles, bbox_to_anchor=(0.5, 0.5), loc="center", fontsize=12)
    fig.savefig(output_path, dpi=dpi)

# %%
def main():
    args = get_args()
    gbk = Genbank(args.input)

    if args.mobsuite_report:
        # Parse the report and keep only contigs classified as 'chromosome'
        df = pd.read_csv(args.mobsuite_report, sep='\t')
        chrom_contigs = df.loc[df["molecule_type"] == "chromosome", "contig_id"].tolist()

        # Validate and extract subset from GBK
        all_features = gbk.get_seqid2features(feature_type=None)
        all_sizes = gbk.get_seqid2size()
        subset_features = {c: all_features[c] for c in chrom_contigs if c in all_features}
        subset_sizes = {c: all_sizes[c] for c in chrom_contigs if c in all_sizes}

        if not subset_features:
            print("Warning: No valid chromosome contigs found in GBK matching the report.")

        output_file = args.output.with_name(f"{args.output.stem}_chromosome.png")
        plot_circos(
            seqid2size=subset_sizes,
            seqid2features=subset_features,
            title=f"{Path(args.input).stem} (chromosome)",
            output_path=output_file,
            mge_elements=args.mge_elements
        )

    else:
        # No report provided: fall back to fragment handling
        all_features = gbk.get_seqid2features(feature_type=None)
        all_sizes = gbk.get_seqid2size()
        filtered = {k: v for k, v in all_sizes.items() if v >= 1000}

        if len(filtered) > 15:
            fig, ax = plt.subplots(figsize=(10, 3))
            ax.text(0.5, 0.5, "Assembly too fragmented. Not plotted.",
                    fontsize=16, ha='center', va='center')
            ax.axis('off')
            fig.savefig(args.output.with_name(f"{args.output.stem}.png"), dpi=300)
            return

        for contig in filtered:
            output_file = args.output.with_name(f"{args.output.stem}_{contig}.png")
            plot_circos(
                seqid2size={contig: all_sizes[contig]},
                seqid2features={contig: all_features[contig]},
                title=f"{Path(args.input).stem} ({contig})",
                output_path=output_file,
                mge_elements=False
            )

# %%
if __name__ == "__main__":
    main()
