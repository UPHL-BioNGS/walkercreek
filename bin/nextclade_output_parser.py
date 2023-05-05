#!/usr/bin/env python
from genericpath import sameopenfile
from importlib.resources import path
from os.path import exists
import argparse
import pandas as pd

# Argument parser: get arguments
parser = argparse.ArgumentParser()
parser.add_argument('--id', required=True, help='The ID of the sample')
args = parser.parse_args()

# Sample id name variable
id_name = args.id

# Read the TSV file from nextclade_run.nf using pandas
tsv_data = pd.read_csv(f"{args.id}", sep='\t')

# If there is only one row, add a row of 'NA' values to ensure that the output files have the same number of rows and the same structure.
if len(tsv_data) == 1:
    tsv_data = tsv_data.append(pd.Series(['NA'] * len(tsv_data.columns)), ignore_index=True)

# Set the index column of the TSV to 'index'
tsv_data = tsv_data.set_index('index')

# Write the output TSV files using pandas
tsv_data.loc[:, ['clade']].to_csv(f"NEXTCLADE_CLADE.tsv", sep='\t', header=False)
tsv_data.loc[:, ['aaSubstitutions']].to_csv(f"NEXTCLADE_AASUBS.tsv", sep='\t', header=False)
tsv_data.loc[:, ['aaDeletions']].to_csv(f"NEXTCLADE_AADELS.tsv", sep='\t', header=False)

# Check if 'lineage' column exists in the input TSV file
if 'lineage' in tsv_data.columns:
    tsv_data.loc[:, ['lineage']].to_csv(f"NEXTCLADE_LINEAGE.tsv", sep='\t', header=False)

# Read in the four output files, handling the case where a file is missing
try:
    df_clade = pd.read_csv(f"NEXTCLADE_CLADE.tsv", sep="\t", header=None, names=["clade"])
except FileNotFoundError:
    print(f"NEXTCLADE_CLADE.tsv not found")
    df_clade = pd.DataFrame(columns=["clade"])
try:
    df_subs = pd.read_csv(f"NEXTCLADE_AASUBS.tsv", sep="\t", header=None, names=["aaSubstitutions"])
except FileNotFoundError:
    print(f"NEXTCLADE_AASUBS.tsv not found")
    df_subs = pd.DataFrame(columns=["aaSubstitutions"])
try:
    df_dels = pd.read_csv(f"NEXTCLADE_AADELS.tsv", sep="\t", header=None, names=["aaDeletions"])
except FileNotFoundError:
    print(f"NEXTCLADE_AADELS.tsv not found")
    df_dels = pd.DataFrame(columns=["aaDeletions"])

# If the 'lineage' column exists in the input TSV file, read in the 'NEXTCLADE_LINEAGE.tsv' output file
if 'lineage' in tsv_data.columns:
    try:
        df_lineage = pd.read_csv("NEXTCLADE_LINEAGE.tsv", sep="\t", header=None, names=["lineage"])
    except FileNotFoundError:
        print("NEXTCLADE_LINEAGE.tsv not found")
        df_lineage = pd.DataFrame(columns=["lineage"])
else:
    df_lineage = pd.DataFrame(columns=["lineage"])

# Concatenate the four dataframes horizontally, excluding the 'lineage' column if it doesn't exist in the input TSV file
df_combined = pd.concat([df_clade, df_subs, df_dels], axis=1)
if 'lineage' in tsv_data.columns:
    df_combined = pd.concat([df_combined, df_lineage], axis=1)

# Write the combined dataframe to a TSV file as final report
df_combined.to_csv(f"NEXTCLADE_REPORT.tsv", sep="\t", index=False)
