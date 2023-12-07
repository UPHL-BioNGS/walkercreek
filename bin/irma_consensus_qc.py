#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqUtils import GC
import sys
import csv


# Function to count specified bases in a given sequence
def count_bases(seq, bases):
    return sum(seq.count(base) for base in bases)


# Function to calculate N50 from IRMA consensus
def calculate_n50(seq_records):
    lengths = sorted([len(rec) for rec in seq_records], reverse=True)
    total_length = sum(lengths)
    count = 0
    for length in lengths:
        count += length
        if count >= total_length / 2:
            return length
    return 0  # In case there are no sequences


# Function to calculate GC content from IRMA consensus
def calculate_gc_content(seq):
    return GC(seq)


# Function to calculate segment counts in IRMA consensus
def calculate_segment_counts(seq_records):
    return len(seq_records)


def main(consensus_fasta, meta_id):
    # Parse the input FASTA file and convert it to a list of records
    consensus_records = list(SeqIO.parse(consensus_fasta, "fasta"))

    # Combine all sequences in the records into a single sequence as assembled gene segments will vary by sample
    combined_seq = "".join([str(record.seq) for record in consensus_records])

    # Calculate counts for different bases
    n_count = count_bases(combined_seq, ["N"])
    actg_count = count_bases(combined_seq, ["A", "C", "T", "G"])
    degenerate_count = len(combined_seq) - actg_count - n_count
    total_count = len(combined_seq)
    segment_count = calculate_segment_counts(consensus_records)
    n50 = calculate_n50(consensus_records)
    gc_content = calculate_gc_content(combined_seq)
    gc_content = round(gc_content, 2)

    # Write the counts to individual output files
    with open("n_count", "w") as f:
        f.write(str(n_count))
    with open("actg_count", "w") as f:
        f.write(str(actg_count))
    with open("degenerate_count", "w") as f:
        f.write(str(degenerate_count))
    with open("total_count", "w") as f:
        f.write(str(total_count))
    with open("segment_count", "w") as f:
        f.write(str(segment_count))
    with open("n50", "w") as f:
        f.write(str(n50))
    with open("gc_content", "w") as f:
        f.write(str(gc_content))

    # Print the counts
    print(f"N count: {n_count}")
    print(f"ACTG count: {actg_count}")
    print(f"Degenerate count: {degenerate_count}")
    print(f"Total count: {total_count}")
    print(f"Segment count: {segment_count}")
    print(f"N50: {n50}")
    print(f"GC content: {gc_content}")

    # Prepare the output file name
    output_tsv = f"{meta_id}.irma_consensus_qc.tsv"

    # Write the computed data to a tab-separated file
    with open(output_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        # Write the header row
        writer.writerow(
            [
                "Sample",
                "IRMA consensus ACTG Count",
                "IRMA consensus Degenerate Count",
                "IRMA consensus N Count",
                "IRMA consensus Total Count",
                "IRMA consensus Segment Count",
                "IRMA consensus N50",
                "IRMA consensus GC Content",
            ]
        )

        # Write the calculated data for the sample
        writer.writerow([meta_id, actg_count, degenerate_count, n_count, total_count, segment_count, n50, gc_content])


# The script's entry point: read command-line arguments and execute the main function
if __name__ == "__main__":
    consensus_fasta = sys.argv[1]
    meta_id = sys.argv[2]
    main(consensus_fasta, meta_id)
