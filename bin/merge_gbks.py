#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from collections import defaultdict
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Merge MGE annotations from GBK files.")
    parser.add_argument("-g", "--genbank_base", required=True, help="Genbank base file with original annotations")
    parser.add_argument("-u", "--updated_gbks", required=False, help="Files with new annotations (list with comma-separated files, no blankspaces allowed)")
    parser.add_argument("-o", "--output", default="merged_output.gbk", help="Outfile name")
    return parser.parse_args()


def feature_key(feature):
    return (feature.type, str(feature.location))


def merge_annotations_multi_record(base_gbk, updated_gbk_list, output_file):
    base_records = SeqIO.to_dict(SeqIO.parse(base_gbk, "genbank"))

    record_features = {rid: {feature_key(f): f for f in rec.features} for rid, rec in base_records.items()}

    for updated_path in updated_gbk_list:
        for rec in SeqIO.parse(updated_path, "genbank"):
            if rec.id not in record_features:
                base_records[rec.id] = rec
                record_features[rec.id] = {
                    feature_key(f): f
                    for f in rec.features
                    if f.qualifiers.get("mge_element", ["no"])[0].lower() == "yes"
                }
            else:
                for feature in rec.features:
                    if feature.qualifiers.get("mge_element", ["no"])[0].lower() == "yes":
                        key = feature_key(feature)
                        record_features[rec.id][key] = feature

    for rec_id, features in record_features.items():
        base_records[rec_id].features = list(features.values())

    SeqIO.write(base_records.values(), output_file, "genbank")
    print(f"[INFO] Final annotation saved in: {output_file}")


if __name__ == "__main__":
    args = parse_args()
    if args.updated_gbks:
        updated_files = args.updated_gbks.split(",")

        for path in [args.genbank_base] + updated_files:
            if not os.path.isfile(path):
                raise FileNotFoundError(f"File not found: {path}")

        merge_annotations_multi_record(args.genbank_base, updated_files, args.output)

    else:
        print("[INFO] No updated GBK files provided, copying base file to output.")
        os.system(f"cp {args.genbank_base} {args.output}")
        print(f"[INFO] Base annotation saved in: {args.output}")
