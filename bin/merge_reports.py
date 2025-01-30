import pandas as pd
import sys


def merge_tsvs(files):
    # Read the first TSV into a DataFrame and ensure 'Sample' is a string
    df = pd.read_csv(files[0], sep="\t")
    if "Sample" in df.columns:
        df["Sample"] = df["Sample"].astype(str)

    # Loop over remaining files and merge
    for file in files[1:]:
        temp_df = pd.read_csv(file, sep="\t")
        if "Sample" in temp_df.columns:
            temp_df["Sample"] = temp_df["Sample"].astype(str)
        df = pd.merge(df, temp_df, on="Sample", how="outer")

    # Round all numeric columns to 2 decimal places
    numeric_cols = df.select_dtypes(include="number").columns
    df[numeric_cols] = df[numeric_cols].round(2)

    return df


def main():
    files = sys.argv[1:]  # List of TSV file paths from the command line arguments
    merged_df = merge_tsvs(files)
    merged_df.to_csv("summary_report.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
