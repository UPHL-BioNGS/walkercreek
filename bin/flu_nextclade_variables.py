#!/usr/bin/env python

import os
import argparse

# Dictionary containing data for various flu subtypes with their associated datasets.
flu_subtypes = {
    "H1N1": {"dataset": "flu_h1n1pdm_ha"},
    "H3N2": {"dataset": "flu_h3n2_ha"},
    "Victoria": {"dataset": "flu_vic_ha"},
    "Yamagata": {"dataset": "flu_yam_ha"},
}


def get_subtype_from_file(file_path, invalid_patterns):
    """
    Reads the subtype from a file if it exists and does not match invalid patterns.
    Returns None if file is absent, empty, or matches an invalid pattern.
    """
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            subtype = f.read().strip()
            # Check for invalid patterns in the subtype
            if subtype and not any(pattern in subtype for pattern in invalid_patterns):
                return subtype
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Outputs the dataset, reference, and tag for the HA gene of a given flu subtype."
    )
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--irma_subtype", required=True, help="Path to IRMA subtype file")
    parser.add_argument("--abricate_subtype", required=True, help="Path to Abricate subtype file")
    args = parser.parse_args()

    # Define patterns that indicate no valid subtype
    invalid_patterns_irma = ["No IRMA subtype", "NoIRMAsubtype"]
    invalid_patterns_abricate = ["No abricate subtype", "Noabricatesubtype"]

    # Attempt to retrieve the flu subtype, prioritizing IRMA if it has a valid subtype
    flu_subtype = get_subtype_from_file(args.irma_subtype, invalid_patterns_irma) or get_subtype_from_file(
        args.abricate_subtype, invalid_patterns_abricate
    )

    # Validate the flu subtype and output the corresponding dataset
    if flu_subtype in flu_subtypes:
        dataset = flu_subtypes[flu_subtype]["dataset"]
        file_path = dataset

        # Write to the file and print information
        with open(file_path, "w") as f:
            f.write(f"{dataset}\n")
            print(f"  {dataset}: {dataset} (output to {file_path})")
    else:
        print(f"Error: No valid flu subtype found for sample '{args.sample}'.")


if __name__ == "__main__":
    main()
