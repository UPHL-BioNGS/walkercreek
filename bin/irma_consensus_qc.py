#!/usr/bin/env python

from Bio import SeqIO
import sys
import csv

def count_bases(seq, bases):
    return sum(seq.count(base) for base in bases)

def calculate_coverage(seq, ref):
    coverage = (len(seq) / len(ref)) * 100
    return round(coverage, 1)

def main(consensus_fasta, reference_fasta, meta_id):
    consensus_records = list(SeqIO.parse(consensus_fasta, "fasta"))
    reference_records = list(SeqIO.parse(reference_fasta, "fasta"))

    for consensus_record in consensus_records:
        consensus_seq = consensus_record.seq

        n_count = count_bases(consensus_seq, ['N'])
        actg_count = count_bases(consensus_seq, ['A', 'C', 'T', 'G'])
        degenerate_count = len(consensus_seq) - actg_count - n_count
        total_count = len(consensus_seq)

        for reference_record in reference_records:
            reference_seq = reference_record.seq

            percent_reference_coverage = calculate_coverage(consensus_seq, reference_seq)

        # Write the counts and coverage values to output files
        with open("n_count", "w") as f:
            f.write(str(n_count))
        with open("actg_count", "w") as f:
            f.write(str(actg_count))
        with open("percent_reference_coverage", "w") as f:
            f.write(str(percent_reference_coverage))
        with open("degenerate_count", "w") as f:
            f.write(str(degenerate_count))
        with open("total_count", "w") as f:
            f.write(str(total_count))

        # Print the counts and coverage information to console output
        print(f"N count: {n_count}")
        print(f"ACTG count: {actg_count}")
        print(f"Degenerate count: {degenerate_count}")
        print(f"Total count: {total_count}")
        print(f"Percent reference coverage: {percent_reference_coverage}%")

    output_tsv = f"{meta_id}.irma_consensus_qc.tsv"
    total_count_sum = 0

    with open(output_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        # Write the header row
        writer.writerow(["Sample Name", "% Reference Coverage Per Gene Segment", "ACTG Count", "Degenerate Count", "N Count", "Total Count"])

        for consensus_record in consensus_records:
            consensus_seq = consensus_record.seq

            n_count = count_bases(consensus_seq, ['N'])
            actg_count = count_bases(consensus_seq, ['A', 'C', 'T', 'G'])
            degenerate_count = len(consensus_seq) - actg_count - n_count
            total_count = len(consensus_seq)
            total_count_sum += total_count

            segment_id = consensus_record.id

            reference_seq = None
            for reference_record in reference_records:
                if reference_record.id == segment_id:
                    reference_seq = reference_record.seq
                    break

            if reference_seq is None:
                print(f"No matching reference segment found for {segment_id}")
                continue

            percent_reference_coverage = calculate_coverage(consensus_seq, reference_seq)

            # Write the row for each record
            writer.writerow([segment_id, percent_reference_coverage, actg_count, degenerate_count, n_count, total_count])

            # Print the counts and coverage information to console output
            print(f"Sample Name: {segment_id}")
            print(f"ACTG count: {actg_count}")
            print(f"Degenerate count: {degenerate_count}")
            print(f"N count: {n_count}")
            print(f"Total count: {total_count}")
            print(f"% Reference Coverage Per Gene Segment: {percent_reference_coverage}%")

    total_percent_ref_coverage = round((total_count_sum / 13500) * 100, 1)

    # Write the total counts and coverage values to output files
    with open(f"{meta_id}.n_count", "w") as f:
        f.write(str(n_count))
    with open(f"{meta_id}.actg_count", "w") as f:
        f.write(str(actg_count))
    with open(f"{meta_id}.percent_reference_coverage", "w") as f:
        f.write(str(percent_reference_coverage))
    with open(f"{meta_id}.degenerate_count", "w") as f:
        f.write(str(degenerate_count))
    with open(f"{meta_id}.total_count", "w") as f:
        f.write(str(total_count))

    with open(output_tsv, "a", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        # Write the total percent reference coverage
        writer.writerow(["Total % Reference Coverage", "", "", "", "", "", total_percent_ref_coverage])

if __name__ == "__main__":
    consensus_fasta = sys.argv[1]
    reference_fasta = sys.argv[2]
    meta_id = sys.argv[3]
    main(consensus_fasta, reference_fasta, meta_id)
