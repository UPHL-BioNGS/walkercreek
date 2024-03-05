import pandas as pd
import sys

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <bam_results> <coverage_results> <output_file>")
        sys.exit(1)
    
    bam_file_path = sys.argv[1]
    coverage_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    # Read the TSV files into DataFrames
    bam_df = pd.read_csv(bam_file_path, sep="\t")
    coverage_df = pd.read_csv(coverage_file_path, sep="\t")

    # Make sure 'Sample' is the first column in both DataFrames
    bam_df = bam_df[['Sample'] + [col for col in bam_df.columns if col != 'Sample']]
    coverage_df = coverage_df[['Sample'] + [col for col in coverage_df.columns if col != 'Sample']]

    # Merge the DataFrames on the 'Sample' column
    merged_df = pd.merge(bam_df, coverage_df, on="Sample", how="outer")

    # Round all numeric columns to 2 decimal places
    merged_df = merged_df.round(2)

    # Make sure 'Sample' is the first column in the merged DataFrame
    merged_df = merged_df[['Sample'] + [col for col in merged_df.columns if col != 'Sample']]

    # Write the merged DataFrame to a TSV file
    merged_df.to_csv(output_file_path, sep="\t", index=False)

if __name__ == "__main__":
    main()
