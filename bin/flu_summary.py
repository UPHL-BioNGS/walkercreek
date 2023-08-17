#!/bin/python

import sys
import pandas as pd

def main(txt_files, tsv_files):
    dataframes = []

    for txt_file in txt_files:
        df = pd.read_csv(txt_file, sep='\t')
        sample_name = get_sample_name(txt_file)
        df['Sample'] = sample_name
        df = format_dataframe(df, sample_name)
        dataframes.append(df)

    for tsv_file in tsv_files:
        df = pd.read_csv(tsv_file, sep='\t')
        sample_name = get_sample_name(tsv_file)
        df['Sample'] = sample_name
        df = format_dataframe(df, sample_name)
        dataframes.append(df)

    if dataframes:
        summary_df = pd.concat(dataframes)
        sorted_df = summary_df.sort_values('Sample')
        sorted_df.to_csv('flu_summary.tsv', sep='\t', index=False)

def format_dataframe(df, sample_name):
    formatted_df = pd.DataFrame(columns=[
        'Sample', 'IRMA type', 'IRMA subtype', 'abricate InsaFlu type', 'abricate InsaFlu subtype',
        'Clade', 'AA Substitutions', 'AA Deletions', 'kraken2 Homo sapiens %', 'kraken2 Influenza A virus %',
        'kraken2 Influenza B virus %',
        'Reads Before Trimming', 'GC Before Trimming', 'Average Q Score Before Trimming',
        'Reads After Trimming', 'Paired Reads After Trimming', 'Unpaired Reads After Trimming',
        'GC After Trimming', 'Average Q Score After Trimming'
    ])

    if 'IRMA_type' in df.columns:
        formatted_df['IRMA type'] = df['IRMA_type']
    if 'IRMA_subtype' in df.columns:
        formatted_df['IRMA subtype'] = df['IRMA_subtype']
    if 'abricate_InsaFlu_type' in df.columns:
        formatted_df['abricate InsaFlu type'] = df['abricate_InsaFlu_type']
    if 'abricate_InsaFlu_subtype' in df.columns:
        formatted_df['abricate InsaFlu subtype'] = df['abricate_InsaFlu_subtype']
    if 'Clade' in df.columns:
        formatted_df['Clade'] = df['Clade']
    if 'aaSubstitutions' in df.columns:
        formatted_df['AA Substitutions'] = df['aaSubstitutions']
    if 'aaDeletions' in df.columns:
        formatted_df['AA Deletions'] = df['aaDeletions']
    if 'Homo sapiens' in df.columns:
        formatted_df['kraken2 Homo sapiens %'] = df['Homo sapiens']
    if 'Influenza A virus' in df.columns:
        formatted_df['kraken2 Influenza A virus %'] = df['Influenza A virus']
    if 'Influenza B virus' in df.columns:
        formatted_df['kraken2 Influenza A virus %'] = df['Influenza B virus']
    if 'Reads Before Trimming' in df.columns:
        formatted_df['Reads Before Trimming'] = df['Reads Before Trimming']
    if 'GC Before Trimming' in df.columns:
        formatted_df['GC Before Trimming'] = df['GC Before Trimming']
    if 'Average Q Score Before Trimming' in df.columns:
        formatted_df['Average Q Score Before Trimming'] = df['Average Q Score Before Trimming']
    if 'Reads After Trimming' in df.columns:
        formatted_df['Reads After Trimming'] = df['Reads After Trimming']
    if 'Paired Reads After Trimming' in df.columns:
        formatted_df['Paired Reads After Trimming'] = df['Paired Reads After Trimming']
    if 'Unpaired Reads After Trimming' in df.columns:
        formatted_df['Unpaired Reads After Trimming'] = df['Unpaired Reads After Trimming']
    if 'GC After Trimming' in df.columns:
        formatted_df['GC After Trimming'] = df['GC After Trimming']
    if 'Average Q Score After Trimming' in df.columns:
        formatted_df['Average Q Score After Trimming'] = df['Average Q Score After Trimming']

    formatted_df['Sample'] = sample_name

    return formatted_df

def get_sample_name(file_path):
    with open(file_path, 'r') as file:
        # Read the first line of the file
        first_line = file.readline()
        # Split the line by tab (\t) to extract the sample name
        columns = first_line.split('\t')
        if len(columns) > 0:
            # Assuming the sample name is in the first column
            sample_name = columns[0].strip()
            return sample_name
    # Return None if the file is empty or cannot be read
    return None

if __name__ == '__main__':
    txt_files = sys.argv[1].split(',') if len(sys.argv) > 1 else []
    tsv_files = sys.argv[2].split(',') if len(sys.argv) > 2 else []
    main(txt_files, tsv_files)
