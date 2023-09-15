#!/usr/bin/env python
import os
from os.path import exists
import argparse
import pandas as pd

# Define the dataset, reference, and tag for each flu subtype
flu_subtypes = {
    "H1N1": {
        "dataset": "flu_h1n1pdm_ha",
        "reference": "CY121680",
        "tag": "2023-04-02T12:00:00Z",
    },
    "H3N2": {
        "dataset": "flu_h3n2_ha",
        "reference": "CY163680",
        "tag": "2023-04-02T12:00:00Z",
    },
    "Victoria": {
        "dataset": "flu_vic_ha",
        "reference": "KX058884",
        "tag": "2023-04-02T12:00:00Z",
    },
    "Yamagata": {
        "dataset": "flu_yam_ha",
        "reference": "JN993010",
        "tag": "2022-07-27T12:00:00Z",
    },
}

def main():
    parser = argparse.ArgumentParser(description="Outputs the dataset, reference, and tag for the HA gene of a given flu subtype.")
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    input_file_path = f"{args.sample}.abricate_flu_subtype.txt"

    if not os.path.exists(input_file_path):
        print(f"Error: Input file '{input_file_path}' does not exist")
        return

    with open(input_file_path, "r") as f:
        flu_subtype = f.read().strip()

    if flu_subtype not in flu_subtypes:
        print(f"Error: Invalid flu subtype '{flu_subtype}' for sample '{args.sample}'")
        return

    for item in ["dataset", "reference", "tag"]:
        file_path = flu_subtypes[flu_subtype][item]
        with open(file_path, "w") as f:
            f.write(f"{flu_subtypes[flu_subtype][item]}\n")
            print(f"  {item}: {flu_subtypes[flu_subtype][item]} (output to {file_path})")

if __name__ == "__main__":
    main()


