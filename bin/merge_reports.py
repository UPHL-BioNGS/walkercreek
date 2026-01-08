import pandas as pd
import sys

def read_tsv(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    # normalize header names
    df.columns = df.columns.str.strip()
    # normalize Sample values if present
    if "Sample" in df.columns:
        df["Sample"] = df["Sample"].astype(str).str.strip()
    return df

def merge_tsvs(files):
    df = read_tsv(files[0])

    for f in files[1:]:
        temp_df = read_tsv(f)
        if "Sample" not in temp_df.columns:
            raise KeyError(f"'Sample' column not found in {f}. Columns: {list(temp_df.columns)}")
        if "Sample" not in df.columns:
            raise KeyError(f"'Sample' column not found in {files[0]}. Columns: {list(df.columns)}")
        df = pd.merge(df, temp_df, on="Sample", how="outer")

    # Round numeric columns
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="ignore")
    numeric_cols = df.select_dtypes(include="number").columns
    df[numeric_cols] = df[numeric_cols].round(2)
    return df

def main():
    files = sys.argv[1:]
    merged_df = merge_tsvs(files)
    merged_df.to_csv("summary_report.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
