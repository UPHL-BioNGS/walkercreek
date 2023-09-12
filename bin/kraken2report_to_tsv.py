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
parser.add_argument("--fluA", help="the target taxonomy for Influenza A virus", default="Influenza A virus")
parser.add_argument("--fluB", help="the target taxonomy for Influenza B virus", default="Influenza B virus")
args = parser.parse_args()

# Sample name variable
sample_name = args.sample

# load the kraken2 report into a pandas dataframe
kraken2_report = pd.read_csv(args.report, sep="\t", header=None)

# filter the dataframe to only rows containing the target taxonomies
human_rows = kraken2_report[kraken2_report[5].str.contains(args.human)]
unclassified_rows = kraken2_report[kraken2_report[5].str.contains(args.unclassified)]
fluA_rows = kraken2_report[kraken2_report[5].str.contains(args.fluA)]
fluB_rows = kraken2_report[kraken2_report[5].str.contains(args.fluB)]

# check if there are any rows for unclassified before extracting the percentage value
if not unclassified_rows.empty:
    temp_percentage_unclassified = unclassified_rows.iloc[0, 0]
    percentage_unclassified = "NA" if temp_percentage_unclassified == 0.00 else temp_percentage_unclassified
else:
    percentage_unclassified = "NA"

# check if there are any rows for Influenza B virus before extracting the percentage value
if not fluB_rows.empty:
    temp_percentage_fluB = fluB_rows.iloc[0, 0]
    percentage_fluB = "NA" if temp_percentage_fluB == 0.00 else temp_percentage_fluB
else:
    percentage_fluB = "NA"

# check if there are any rows for Influenza A virus before extracting the percentage value
if not fluA_rows.empty:
    temp_percentage_fluA = fluA_rows.iloc[0, 0]
    percentage_fluA = "NA" if temp_percentage_fluA == 0.00 else temp_percentage_fluA
else:
    percentage_fluA = "NA"

# check if there are any rows for Homo sapiens before extracting the percentage value
if not human_rows.empty:
    temp_percentage_human = human_rows.iloc[0, 0]
    percentage_human = "NA" if temp_percentage_human == 0.00 else temp_percentage_human
else:
    percentage_human = "NA"

# write the output to the output file
with open(f"{sample_name}_read_percentages.txt", "w") as outfile:
    outfile.write(f"{args.sample}\t{percentage_human}\t{percentage_fluA}\t{percentage_fluB}\t{percentage_unclassified}")

