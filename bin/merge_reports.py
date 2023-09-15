import pandas as pd
import sys

def merge_tsvs(files):
    # Read the first TSV into a DataFrame
    df = pd.read_csv(files[0], sep="\t")

    # Loop over remaining files and merge
    for file in files[1:]:
        temp_df = pd.read_csv(file, sep="\t")
        df = pd.merge(df, temp_df, on='Sample', how='outer')

    return df

def main():
    files = sys.argv[1:]  # List of TSV file paths from the command line arguments
    merged_df = merge_tsvs(files)
    merged_df.to_csv('summary_report.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()

