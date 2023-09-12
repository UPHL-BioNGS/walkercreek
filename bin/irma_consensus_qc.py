#!/usr/bin/env python

from Bio import SeqIO
import sys
import csv

# Function to count specified bases in a given sequence
def count_bases(seq, bases):
    return sum(seq.count(base) for base in bases)

# Function to calculate coverage of the genome based on the sequence length
def calculate_coverage(seq, genome_length):
    coverage = (len(seq) / genome_length) * 100
    return round(coverage, 1)

def main(consensus_fasta, meta_id):
    # Parse the input FASTA file and convert it to a list of records
    consensus_records = list(SeqIO.parse(consensus_fasta, "fasta"))

    # Combine all sequences in the records into a single sequence as assembled gene segments will vary by sample
    combined_seq = ''.join([str(record.seq) for record in consensus_records])

    # Calculate counts for different bases
    n_count = count_bases(combined_seq, ['N'])
    actg_count = count_bases(combined_seq, ['A', 'C', 'T', 'G'])
    degenerate_count = len(combined_seq) - actg_count - n_count
    total_count = len(combined_seq)

    # Calculate percent coverage using a predefined genome length due to IRMA's
    percent_coverage = calculate_coverage(combined_seq, 14200) #14200 used as a avg combo for fluA and fluB genome lengths

    # Write the counts and coverage values to individual output files
    with open("n_count", "w") as f:
        f.write(str(n_count))
    with open("actg_count", "w") as f:
        f.write(str(actg_count))
    with open("percent_coverage", "w") as f:
        f.write(str(percent_coverage))
    with open("degenerate_count", "w") as f:
        f.write(str(degenerate_count))
    with open("total_count", "w") as f:
        f.write(str(total_count))

    # Print the counts and coverage information to console
    print(f"N count: {n_count}")
    print(f"ACTG count: {actg_count}")
    print(f"Degenerate count: {degenerate_count}")
    print(f"Total count: {total_count}")
    print(f"Percent coverage: {percent_coverage}%")

    # Prepare the output file name
    output_tsv = f"{meta_id}.irma_consensus_qc.tsv"

    # Write the computed data to a tab-separated file
    with open(output_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        # Write the header row
        writer.writerow(["Sample Name", "% Coverage", "ACTG Count", "Degenerate Count", "N Count", "Total Count"])

        # Write the calculated data for the sample
        writer.writerow([meta_id, percent_coverage, actg_count, degenerate_count, n_count, total_count])

# The script's entry point: read command-line arguments and execute the main function
if __name__ == "__main__":
    consensus_fasta = sys.argv[1]
    meta_id = sys.argv[2]
    main(consensus_fasta, meta_id)
