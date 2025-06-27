#!/usr/bin/env python3
import argparse
import logging
from pathlib import Path
from typing import List, Set

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

# -----------------------------------------------------------------------------
# CONFIGURE LOGGING
# -----------------------------------------------------------------------------
logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# ARG PARSING
# -----------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Merge MGE summary TSVs into one table"
    )
    p.add_argument("-s", "--sample", required=True,
                   help="Sample name (goes in the Sample column)")
    p.add_argument("-g", "--genbank", required=True, type=Path,
                   help="Assembly GenBank file (.gbk)")
    p.add_argument("-t", "--tables",
                   help="Comma-separated list of *_summary.tsv files")
    p.add_argument("-o", "--out", required=True, type=Path,
                   help="Output prefix (writes <out>.tsv)")
    return p.parse_args()


# -----------------------------------------------------------------------------
# UTILS
# -----------------------------------------------------------------------------
def load_contigs(gbk_path: Path) -> Set[str]:
    """Extract contig IDs from a GenBank file."""
    contigs: Set[str] = set()
    for rec in SeqIO.parse(str(gbk_path), "genbank"):
        contigs.add(rec.id)
    logger.info(f"Found {len(contigs)} contigs in {gbk_path.name}")
    return contigs


def infer_mge_type(path: Path) -> str:
    """Infer MGE type from filename, e.g. 'prophage' from 'prophage_summary.tsv'."""
    name = path.stem  # e.g. "prophage_summary"
    if name.endswith("_summary"):
        return name[: -len("_summary")]
    return name


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
def main():
    args = parse_args()

    # 1) load contigs for validation
    try:
        gbk_records = list(SeqIO.parse(str(args.genbank), "genbank"))
        valid_contigs = {rec.id for rec in gbk_records}
    except Exception as e:
        logger.critical(f"Could not read GenBank {args.genbank}: {e}")
        raise

    # 2) process each summary table
    merged: List[pd.DataFrame] = []
    mge_features_by_contig = {rec.id: [] for rec in gbk_records}

    if args.tables is None:
        logger.warning("No tables provided.")
        mge_features = []
    else:
        for tbl_path in map(Path, args.tables.split(",")):
            mge = infer_mge_type(tbl_path)
            logger.info(f"Processing {tbl_path.name} as MGE='{mge}'")

            # read with pandas
            df = pd.read_table(tbl_path, dtype=str)

            # warn if contigs not in GenBank
            if "Contig" in df.columns:
                missing = set(df["Contig"]) - valid_contigs
                if missing:
                    logger.warning(
                        f"{len(missing)} contigs in {tbl_path.name} not found in GenBank: {', '.join(list(missing))}"
                    )

            # pick and rename the needed columns
            df2 = df[["Contig", "Start", "End", "Length", "Name", "AMR", "VF"]].copy()
            df2["Sample"] = args.sample
            df2["MGE"] = mge
            for _, row in df2.iterrows():
                try:
                    contig = row["Contig"]
                    start = int(row["Start"])
                    end = int(row["End"])
                    name = row["Name"]

                    feature = SeqFeature(
                        FeatureLocation(start - 1, end),  # GenBank is 0-based start
                        type="MGE",
                        qualifiers={
                            "type": mge,
                            "name": name
                        }
                    )
                    mge_features_by_contig[contig].append(feature)
                except Exception as e:
                    logger.error(f"Could not create MGE feature from row: {row.to_dict()}\n{e}")

            # reorder
            df2 = df2[["Sample", "Contig", "Start", "End", "Length", "MGE", "Name", "AMR", "VF"]]
            merged.append(df2)

        # 3) concatenate and write out
        out_df = pd.concat(merged, ignore_index=True)
        out_df.sort_values(
                            by=["Contig", "Start"],
                            ascending=[True, True],
                            inplace=True
                            )
        out_file = args.out.with_suffix(".tsv")
        out_df.to_csv(out_file, sep="\t", index=False)
        logger.info(f"Written merged table to {out_file}")


    # 4) write GenBank with added MGE features
    updated_records = []
    for rec in gbk_records:
        orig_features = rec.features
        mge_features = sorted(mge_features_by_contig.get(rec.id, []), key=lambda f: f.location.start)

        # merge features keeping the order by start position
        combined = sorted(orig_features + mge_features, key=lambda f: f.location.start)
        rec.features = combined
        updated_records.append(rec)

    gbk_out = args.out.with_suffix(".gbk")
    with gbk_out.open("w") as f:
        SeqIO.write(updated_records, f, "genbank")
    logger.info(f"Written updated GenBank with MGE features to {gbk_out}")

if __name__ == "__main__":
    main()
