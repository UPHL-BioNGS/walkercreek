#!/usr/bin/env python
from genericpath import sameopenfile
from importlib.resources import path
from os.path import exists
import argparse
import pandas as pd
import re

# Argument parser: get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--id", required=True, help="The ID of the sample")
args = parser.parse_args()

# Sample id name variable
id_name = args.id

# Read the TSV file from nextclade_run.nf using pandas
tsv_data = pd.read_csv(f"{args.id}", sep="\t")

# Set the index column of the TSV to 'index'
tsv_data = tsv_data.set_index("index")

# Convert coverage to percentage
tsv_data["coverage"] = (tsv_data["coverage"] * 100).round(2)

# Write the output TSV files using pandas
tsv_data.loc[:, ["clade"]].to_csv(f"NEXTCLADE_CLADE.tsv", sep="\t", header=False)

# Check if 'short_clade' and 'subclade' columns exist in the input TSV file
if "short_clade" in tsv_data.columns:
    df_short_clade = tsv_data.loc[:, ["short_clade"]]
else:
    df_short_clade = pd.DataFrame(columns=["short_clade"])

if "subclade" in tsv_data.columns:
    df_subclade = tsv_data.loc[:, ["subclade"]]
else:
    df_subclade = pd.DataFrame(columns=["subclade"])

# If the 'lineage' column exists in the input TSV file, read in the 'NEXTCLADE_LINEAGE.tsv' output file
if "lineage" in tsv_data.columns:
    try:
        df_lineage = pd.read_csv("NEXTCLADE_LINEAGE.tsv", sep="\t", header=None, names=["lineage"])
    except FileNotFoundError:
        print("NEXTCLADE_LINEAGE.tsv not found")
        df_lineage = pd.DataFrame(columns=["lineage"])
else:
    df_lineage = pd.DataFrame(columns=["lineage"])

# Read in the clade output file, handling the case where a file is missing
try:
    df_clade = pd.read_csv(f"NEXTCLADE_CLADE.tsv", sep="\t", header=None, names=["clade"])
except FileNotFoundError:
    print(f"NEXTCLADE_CLADE.tsv not found")
    df_clade = pd.DataFrame(columns=["clade"])

# Extract the desired columns from the original TSV data
df_qc_score = tsv_data.loc[:, ["qc.overallScore"]]
df_qc_status = tsv_data.loc[:, ["qc.overallStatus"]]
df_substitutions = tsv_data.loc[:, ["totalSubstitutions"]]
df_coverage = tsv_data.loc[:, ["coverage"]]

# Remove file extension from the sample name
sample_name = id_name.rsplit(".", 1)[0]

# Create a DataFrame with the sample name
df_sample = pd.DataFrame({"Sample": [sample_name]})

# Concatenate all the dataframes horizontally
df_combined = pd.concat(
    [df_sample, df_clade, df_short_clade, df_subclade, df_qc_score, df_qc_status, df_substitutions, df_coverage], axis=1
)

if "lineage" in tsv_data.columns:
    df_combined = pd.concat([df_combined, df_lineage], axis=1)

# Create a dictionary to rename columns with the desired prefix
column_prefix = "Nextclade_"
column_mapping = {col: f"{column_prefix}{col}" for col in df_combined.columns}

# Remove the prefix for the df_sample, df_clade, and df_subclade columns
column_mapping["Sample"] = "Sample"
column_mapping["clade"] = "clade"
column_mapping["subclade"] = "subclade"
column_mapping["short_clade"] = "short_clade"

# Rename the columns using the dictionary
df_combined = df_combined.rename(columns=column_mapping)

# Write the combined dataframe to a TSV file as final report
df_combined.to_csv(f"{sample_name}.nextclade_report.tsv", sep="\t", index=False)
