#!/usr/bin/env python3

import os
import sys
import logging
import argparse
import math
import pandas as pd
import numpy as np
from pathlib import Path


def parse_gff(file_path):
    cds_dict = {}

    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])

            attributes = fields[8].split(";")
            attr_dict = {}
            for attr in attributes:
                attr = attr.strip()
                if attr:
                    key, value = attr.split("=")
                    attr_dict[key] = value

            if feature_type == "CDS":
                product = attr_dict.get("product", "")
                cds_dict[(start, end)] = product

    return cds_dict


def find_matching_cds_range(locus, cds_dict):
    matching_cds_ranges = []
    for key in cds_dict.keys():
        if isinstance(key, tuple) and len(key) == 2:
            if locus >= key[0] and locus <= key[1]:
                matching_cds_ranges.append(key)

    if len(matching_cds_ranges) > 0:
        # return the first matching cds range, there should only be one except for some regions in M2
        return matching_cds_ranges[0]
    else:
        return np.nan


def find_gene_product(locus, cds_dict):
    cds_range = find_matching_cds_range(locus, cds_dict)
    if isinstance(cds_range, tuple):
        return cds_dict[cds_range]
    else:
        if math.isnan(cds_range):
            return np.nan


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Edit iVar variants tsv output to specific format",
        epilog="Example: python edit_ivar_variants.py ivar_variants.tsv ref.gff ivar_variants.edited",
    )
    parser.add_argument("ivar_variants", metavar="IVAR_VARIANTS", help="iVar variants tsv", type=Path)
    parser.add_argument("gff", metavar="GFF", help="Reference gff file", type=Path)
    parser.add_argument("file_out_name", metavar="FILE_OUT_NAME", help="Output file name", type=str)
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.ivar_variants.is_file():
        logging.error(f"The given input file {args.ivar_variants} was not found!")
        sys.exit(2)
    if not args.gff.is_file():
        logging.error(f"The given input file {args.gff} was not found!")
        sys.exit(2)

    sample_name = args.ivar_variants.name.split(".")[0]

    cds_dict = parse_gff(args.gff)
    df_variants = pd.read_csv(args.ivar_variants, sep="\t")
    df_variants["GENE_PRODUCT"] = df_variants["POS"].apply(lambda x: find_gene_product(x, cds_dict))
    df_variants["POS_AA"] = df_variants["POS_AA"].astype("Int64")
    df_variants["SNPID"] = df_variants["REF"] + df_variants["POS"].astype(int).astype(str) + df_variants["ALT"]

    df_variants_new = df_variants[
        ["GENE_PRODUCT", "POS_AA", "REF_AA", "ALT_AA", "TOTAL_DP", "ALT_DP", "ALT_FREQ", "POS", "SNPID"]
    ]
    df_variants_new.insert(0, "SAMPLE", sample_name)
    df_variants_new.to_csv(f"{args.file_out_name}.tsv", sep="\t", index=False)

    # df_variants_nonsynonamous= \
    #    df_variants_new[(df_variants_new['REF_AA'] != df_variants_new['ALT_AA']) \
    #    & df_variants_new['REF_AA'].notnull() \
    #    & df_variants_new['ALT_AA'].notnull()]
    # df_variants_nonsynonamous.to_csv(f"{sample_name}.variants.nonsynonamous.tsv",sep='\t',index=False)


if __name__ == "__main__":
    sys.exit(main())
