#!/usr/bin/env python3

"""
Script to filter a FASTA .faa file by contigs specified from a GenBank file.
"""

import argparse
from pathlib import Path
from typing import Set, Iterator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging

def get_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        prog="faa_contig_filter.py",
        description="Filter proteins in an FAA file by contigs."
    )

    parser.add_argument(
        "-f", "--faa", required=True,
        type=Path,
        help="Required. Input FAA file (fasta)"
    )
    parser.add_argument(
        "-g", "--gbk", required=True,
        type=Path,
        help="Required. Input GenBank file"
    )
    parser.add_argument(
        "-c", "--contigs", required=True,
        type=Path,
        help="Required. File listing contig names, one per line"
    )
    parser.add_argument(
        "-o", "--out", required=True,
        type=Path,
        help="Required. Output directory"
    )

    return parser.parse_args()

def read_contigs(path: Path) -> Set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}

def get_locus_tags(gbk_path: Path, contigs: Set[str]) -> Set[str]:
    """
    Parse GenBank, collect locus_tags for specified contigs
    """
    tags = set()
    for record in SeqIO.parse(str(gbk_path), "genbank"):
        if record.name in contigs:
            for feat in record.features:
                if feat.type == "gene" and "locus_tag" in feat.qualifiers:
                    tags.add(feat.qualifiers["locus_tag"][0])
    return tags

def filter_faa(faa_path: Path, locus_tags: Set[str]) -> Iterator[SeqRecord]:
    """
    Yield SeqRecord entries whose ID matches any of the locus_tags
    """
    for rec in SeqIO.parse(str(faa_path), "fasta"):
        if rec.id.split()[0] in locus_tags:
            yield rec

def main():
    args = get_args()

    # Prepare
    args.out.mkdir(parents=True, exist_ok=True)
    out_file = args.out / f"{args.faa.stem}_filtered.faa"
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    # Read contigs
    contig_set = read_contigs(args.contigs)
    logging.info(f"Loaded {len(contig_set)} contig IDs")

    # Get locus tags
    locus_tags = get_locus_tags(args.gbk, contig_set)
    logging.info(f"Found {len(locus_tags)} locus tags in GenBank file")

    # Filter and write
    filtered = list(filter_faa(args.faa, locus_tags))
    SeqIO.write(filtered, str(out_file), "fasta")
    logging.info(f"Wrote {len(filtered)} records to: {out_file}")

if __name__ == "__main__":
    main()
