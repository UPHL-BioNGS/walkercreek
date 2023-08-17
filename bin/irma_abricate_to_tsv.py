#!/usr/bin/env python

import argparse
import pandas as pd
import os

# Argument parser: get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--sample")
parser.add_argument("--irma_type", type=argparse.FileType("r"), default=None)
parser.add_argument("--irma_subtype", type=argparse.FileType("r"), default=None)
parser.add_argument("--abricate_type", type=argparse.FileType("r"), default=None)
parser.add_argument("--abricate_subtype", type=argparse.FileType("r"), default=None)
args = parser.parse_args()

# Sample name variable
sample = args.sample.strip() if args.sample else ""  # Remove any leading/trailing whitespaces

# Read the first line from each file and remove newline characters
irma_type = args.irma_type.readline().strip() if args.irma_type else ""
irma_subtype = args.irma_subtype.readline().strip() if args.irma_subtype else ""
abricate_type = args.abricate_type.readline().strip() if args.abricate_type else ""
abricate_subtype = args.abricate_subtype.readline().strip() if args.abricate_subtype else ""

# Preparing output list with variables and then reformat into a string
output_list = [
    sample,
    irma_type,
    irma_subtype,
    abricate_type,
    abricate_subtype
]

# Create tab-delimited string for report generation
print('\t'.join(output_list))
