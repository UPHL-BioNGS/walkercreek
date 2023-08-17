#!/usr/bin/env python

from genericpath import sameopenfile
from importlib.resources import path
import argparse
import pandas as pd

# Argument parser: get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--sample")
parser.add_argument("--clade", type=str, help='Path to the clade file')
parser.add_argument("--nextclade_qc_score", type=str, help='Path to the nextclade QC score file')
parser.add_argument("--nextclade_qc_status", type=str, help='Path to the nextclade QC status file')
parser.add_argument("--nextclade_total_substitutions", type=str, help='Path to the nextclade total substitutions file')
parser.add_argument("--nextclade_gene_segment_coverage", type=str, help='Path to the nextclade gene segment coverage file')
parser.add_argument("--nextclade_substitutions", type=str, help='Path to the nextclade substitutions file')
parser.add_argument("--tsv_file", type=str, help='Path to the TSV file')
args = parser.parse_args()

# Sample name variable
sample = args.sample

# Read the TSV file using pandas or set it as an empty DataFrame
tsv_data = pd.DataFrame()
if args.tsv_file:
    tsv_data = pd.read_csv(args.tsv_file, sep='\t')

# Access the values from the TSV file using the column names or set them as empty strings
clade_value = ""
nextclade_qc_score_value = ""
nextclade_qc_status_value = ""
nextclade_total_substitutions_value = ""
nextclade_gene_segment_coverage_value = ""
nextclade_substitutions_value = ""
if not tsv_data.empty:
    if args.clade and args.clade in tsv_data.columns:
        clade_value = tsv_data[args.clade].iloc[0, 1]
    if args.nextclade_qc_score and args.nextclade_qc_score in tsv_data.columns:
        nextclade_qc_score_value = tsv_data[args.nextclade_qc_score].iloc[0, 4]
    if args.nextclade_qc_status and args.nextclade_qc_status in tsv_data.columns:
        nextclade_qc_status_value = tsv_data[args.nextclade_qc_status].iloc[0, 5]
    if args.nextclade_total_substitutions and args.nextclade_total_substitutions in tsv_data.columns:
        nextclade_total_substitutions_value = tsv_data[args.nextclade_total_substitutions].iloc[0, 6]
    if args.nextclade_gene_segment_coverage and args.nextclade_gene_segment_coverage in tsv_data.columns:
        nextclade_gene_segment_coverage_value = tsv_data[args.nextclade_gene_segment_coverage].iloc[0, 7]
    if args.nextclade_substitutions and args.nextclade_substitutions in tsv_data.columns:
        nextclade_substitutions_value = tsv_data[args.nextclade_substitutions].iloc[0, 2]

# Preparing output list with variables and then reformatting into a string
output_list = [
    sample,
    str(clade_value),
    str(nextclade_qc_score_value),
    str(nextclade_qc_status_value),
    str(nextclade_total_substitutions_value),
    str(nextclade_gene_segment_coverage_value),
    str(nextclade_substitutions_value),
]

# Create tab-delimited string for report generation
print('\t'.join(output_list))



