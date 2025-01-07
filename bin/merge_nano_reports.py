import pandas as pd
import sys


def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <raw_results> <filtered_results> <output_file>")
        sys.exit(1)

    raw_file_path = sys.argv[1]
    filtered_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    # Read the TSV files into DataFrames
    raw_df = pd.read_csv(raw_file_path, sep="\t")
    filtered_df = pd.read_csv(filtered_file_path, sep="\t")

    # Make sure 'Sample' is the first column in both DataFrames
    raw_df = raw_df[["Sample"] + [col for col in raw_df.columns if col != "Sample"]]
    filtered_df = filtered_df[["Sample"] + [col for col in filtered_df.columns if col != "Sample"]]

    # Merge the DataFrames on the 'Sample' column
    merged_df = pd.merge(raw_df, filtered_df, on="Sample", how="outer")

    # Round all numeric columns to 2 decimal places
    merged_df = merged_df.round(2)

    # Make sure 'Sample' is the first column in the merged DataFrame
    merged_df = merged_df[["Sample"] + [col for col in merged_df.columns if col != "Sample"]]

    # Write the merged DataFrame to a TSV file
    merged_df.to_csv(output_file_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
