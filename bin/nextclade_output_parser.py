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

# If the 'lineage' column exists in the input TSV file, read in the 'NEXTCLADE_LINEAGE.tsv' output file
if 'lineage' in tsv_data.columns:
    try:
        df_lineage = pd.read_csv("NEXTCLADE_LINEAGE.tsv", sep="\t", header=None, names=["lineage"])
    except FileNotFoundError:
        print("NEXTCLADE_LINEAGE.tsv not found")
        df_lineage = pd.DataFrame(columns=["lineage"])
else:
    df_lineage = pd.DataFrame(columns=["lineage"])

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

# Extract the desired columns from the original TSV data
df_qc_score = tsv_data.loc[:, ['qc.overallScore']]
df_qc_status = tsv_data.loc[:, ['qc.overallStatus']]
df_substitutions = tsv_data.loc[:, ['totalSubstitutions']]
df_coverage = tsv_data.loc[:, ['coverage']]

# Remove file extension from the sample name
sample_name = id_name.rsplit('.', 1)[0]

# Create a DataFrame with the sample name
df_sample = pd.DataFrame({'Sample': [sample_name]})

# Concatenate all the dataframes horizontally
df_combined = pd.concat([df_sample, df_clade, df_subs, df_dels, df_qc_score, df_qc_status, df_substitutions, df_coverage], axis=1)

if 'lineage' in tsv_data.columns:
    df_combined = pd.concat([df_combined, df_lineage], axis=1)

# Write the combined dataframe to a TSV file as final report
df_combined.to_csv(f"{sample_name}.nextclade_report.tsv", sep="\t", index=False)