#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser(description="Parse NanoPlot stats files and produce a summary line.")
parser.add_argument("--sample", required=True, help="Sample name")
parser.add_argument("--filtered", type=argparse.FileType("r"))

args = parser.parse_args()

sample = args.sample


# Function to parse NanoStats file and extract metrics
def parse_nanostats(file_handle):
    mean_length = None
    mean_quality = None
    num_reads = None
    total_bases = None

    for line in file_handle:
        line = line.strip()
        if line.startswith("Mean read length:"):
            mean_length = line.split(":")[1].strip()
        elif line.startswith("Mean read quality:"):
            mean_quality = line.split(":")[1].strip()
        elif line.startswith("Number of reads:"):
            num_reads = line.split(":")[1].strip()
        elif line.startswith("Total bases:"):
            total_bases = line.split(":")[1].strip()

    return mean_length, mean_quality, num_reads, total_bases


filt_mean_length, filt_mean_quality, filt_num_reads, filt_total_bases = parse_nanostats(args.filtered)

# Create a single line of output:
# Sample, Raw stats (4 fields), Post-filtering stats (4 fields), Filtered stats (4 fields)
output_list = [sample, str(filt_mean_length), str(filt_mean_quality), str(filt_num_reads), str(filt_total_bases)]

print("\t".join(output_list))
