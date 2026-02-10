#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV (wide segment metrics)")
    ap.add_argument("--out", required=True, help="Output TSV (filtered)")
    ap.add_argument(
        "--segments",
        default="A_HA,A_NA,B_HA,B_NA",
        help="Comma-separated segments to keep (default: A_HA,A_NA,B_HA,B_NA)",
    )
    ap.add_argument(
        "--keep",
        default="mapped_reads,mean_depth,percent_coverage,reference_length,seq_length",
        help="Comma-separated metric suffixes to keep (default: mapped_reads,mean_depth,percent_coverage,reference_length,seq_length)",
    )
    ap.add_argument(
        "--require-sample",
        action="store_true",
        help="Fail if 'Sample' column is missing (default: no)",
    )
    args = ap.parse_args()

    segs = [s.strip() for s in args.segments.split(",") if s.strip()]
    metrics = [m.strip() for m in args.keep.split(",") if m.strip()]

    df = pd.read_csv(args.inp, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()

    if "Sample" not in df.columns:
        if args.require_sample:
            print("ERROR: Input TSV missing 'Sample' column", file=sys.stderr)
            sys.exit(1)
        else:
            # Still write something sane
            df.to_csv(args.out, sep="\t", index=False)
            return

    # Build desired columns in stable order
    desired = ["Sample"]
    for seg in segs:
        for m in metrics:
            col = f"{seg}_{m}"
            if col in df.columns:
                desired.append(col)

    out = df.loc[:, desired].copy()
    out.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
