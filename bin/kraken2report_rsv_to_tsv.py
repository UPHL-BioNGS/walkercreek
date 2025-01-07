#!/usr/bin/env python
from genericpath import sameopenfile
from importlib.resources import path
from os.path import exists
import argparse
import glob
import pandas as pd

# Argument parser: get arguments from Kraken2 report text file
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sample")
parser.add_argument("--report", type=argparse.FileType("r"))
parser.add_argument("--human", help="the target taxonomy for human", default="Homo sapiens")
parser.add_argument("--unclassified", help="the target taxonomy for unclassified", default="unclassified")
parser.add_argument("--rsv", help="the target taxonomy for RSV virus", default="Respiratory syncytial virus")
args = parser.parse_args()

# Sample name variable
sample_name = args.sample

# load the kraken2 report into a pandas dataframe
kraken2_report = pd.read_csv(args.report, sep="\t", header=None)

# filter the dataframe to only rows containing the target taxonomies
human_rows = kraken2_report[kraken2_report[5].str.contains(args.human)]
unclassified_rows = kraken2_report[kraken2_report[5].str.contains(args.unclassified)]
rsv_rows = kraken2_report[kraken2_report[5].str.contains(args.rsv)]

# check if there are any rows for unclassified before extracting the percentage value
if not unclassified_rows.empty:
    temp_percentage_unclassified = unclassified_rows.iloc[0, 0]
    percentage_unclassified = "NA" if temp_percentage_unclassified == 0.00 else temp_percentage_unclassified
else:
    percentage_unclassified = "NA"

# check if there are any rows for RSV A virus before extracting the percentage value
if not rsv_rows.empty:
    temp_percentage_rsv = rsv_rows.iloc[0, 0]
    percentage_rsv = "NA" if temp_percentage_rsv == 0.00 else temp_percentage_rsv
else:
    percentage_rsv = "NA"

# check if there are any rows for Homo sapiens before extracting the percentage value
if not human_rows.empty:
    temp_percentage_human = human_rows.iloc[0, 0]
    percentage_human = "NA" if temp_percentage_human == 0.00 else temp_percentage_human
else:
    percentage_human = "NA"

# write the output to the output file
with open(f"{sample_name}_read_percentages.txt", "w") as outfile:
    outfile.write(f"{args.sample}\t{percentage_human}\t{percentage_rsv}\t{percentage_unclassified}")
