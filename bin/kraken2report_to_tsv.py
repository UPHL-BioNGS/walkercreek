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
parser.add_argument("--fluA", help="the target taxonomy for Influenza A virus", default="Influenza A virus")
parser.add_argument("--fluB", help="the target taxonomy for Influenza B virus", default="Influenza B virus")
args = parser.parse_args()

# Sample name variable
sample_name = args.sample

# load the kraken2 report into a pandas dataframe
kraken2_report = pd.read_csv(args.report, sep="\t", header=None)

# filter the dataframe to only rows containing the target taxonomies
human_rows = kraken2_report[kraken2_report[5].str.contains(args.human)]
fluA_rows = kraken2_report[kraken2_report[5].str.contains(args.fluA)]
fluB_rows = kraken2_report[kraken2_report[5].str.contains(args.fluB)]

# check if there are any rows for Influenza B virus before extracting the percentage value
if not fluB_rows.empty:
    percentage_fluB = fluB_rows.iloc[0,0]
else:
    percentage_fluB = "NA"


# extract the percentage values from the first column of the filtered rows
percentage_human = human_rows.iloc[0, 0] if not human_rows.empty else "NA"
percentage_fluA = fluA_rows.iloc[0, 0] if not fluA_rows.empty else "NA"

# create a tab-delimited string for report with headers
output = f"File\t{args.human}\t{args.fluA}\t{args.fluB}\n"
output += f"{args.report.name}\t{percentage_human}\t{percentage_fluA}\t{percentage_fluB}"

# write the output to the output file
with open(f"{sample_name}_read_percentages.tsv", "w") as outfile:
    outfile.write(output)

# Generate the summary file
summary_files = glob.glob("*_read_percentages.tsv")
summary_data = []

for file in summary_files:
    sample_id = file.replace("_read_percentages.tsv", "")
    with open(file, "r") as infile:
        percentages = infile.read().strip().split("\n")[1]
        summary_data.append(f"{sample_id}\t{percentages}")

summary_data.sort()

with open("kraken2_report_summary.tsv", "w") as summary_file:
    summary_file.write("Sample\tHomo sapiens\tInfluenza A virus\tInfluenza B virus\n")
    summary_file.write("\n".join(summary_data))
