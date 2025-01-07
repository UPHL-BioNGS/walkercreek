#!/usr/bin/env python

# Import required modules
import os
from os.path import exists
import argparse

# Dictionary containing data for various rsv subtypes with their associated datasets.
rsv_subtypes = {
    "AD": {"dataset": "rsv_a"},
    "BD": {"dataset": "rsv_b"},
    "A": {"dataset": "rsv_a"},
    "B": {"dataset": "rsv_b"},
}


def main():
    # Set up an argument parser to accept input arguments for the script
    parser = argparse.ArgumentParser(description="Outputs the dataset, reference, and tag for a given rsv subtype.")
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    # Construct the path of the input file based on the provided sample name
    input_file_path = f"{args.sample}.irma_subtype.txt"

    # Check if the input file exists
    if not os.path.exists(input_file_path):
        print(f"Error: Input file '{input_file_path}' does not exist")
        return

    # Read the rsv subtype from the input file
    with open(input_file_path, "r") as f:
        rsv_subtype = f.read().strip()

    # Verify that the read subtype is valid and present in the rsv_subtypes dictionary
    if rsv_subtype not in rsv_subtypes:
        print(f"Error: Invalid rsv subtype '{rsv_subtype}' for sample '{args.sample}'")
        return

    # Prepare the dataset and file_path
    dataset = rsv_subtypes[rsv_subtype]["dataset"]
    file_path = dataset

    # Write to the file and print information
    with open(file_path, "w") as f:
        f.write(f"{dataset}\n")
        print(f"  {dataset}: {dataset} (output to {file_path})")


if __name__ == "__main__":
    main()
