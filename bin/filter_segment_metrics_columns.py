#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

DEFAULT_FLU = "A_HA,A_NA,B_HA,B_NA"
DEFAULT_RSV = "RSV_A,RSV_AD,RSV_B,RSV_BD"


def detect_platform_from_columns(cols):
    """
    Detect platform by looking at existing wide-format columns.
    """
    cols = [c.strip() for c in cols]
    if any(c.startswith("RSV_") for c in cols):
        return "rsv"
    if any(c.startswith("A_") or c.startswith("B_") for c in cols):
        return "flu"
    return "flu"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV (wide segment metrics)")
    ap.add_argument("--out", required=True, help="Output TSV (filtered)")
    ap.add_argument(
        "--segments",
        default="",
        help="Comma-separated segments to keep. If empty, chosen by --platform or auto-detect.",
    )
    ap.add_argument(
        "--keep",
        default="mapped_reads,mean_depth,percent_coverage,reference_length,seq_length",
        help="Comma-separated metric suffixes to keep",
    )
    ap.add_argument(
        "--platform",
        default="auto",
        choices=["auto", "flu", "rsv"],
        help="Controls default segment set when --segments not provided. Default: auto.",
    )
    ap.add_argument(
        "--require-sample",
        action="store_true",
        help="Fail if 'Sample' column is missing (default: no)",
    )
    args = ap.parse_args()

    df = pd.read_csv(args.inp, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()

    if "Sample" not in df.columns:
        if args.require_sample:
            print("ERROR: Input TSV missing 'Sample' column", file=sys.stderr)
            sys.exit(1)
        df.to_csv(args.out, sep="\t", index=False)
        return

    # Decide platform if needed
    platform = args.platform
    if platform == "auto":
        platform = detect_platform_from_columns(df.columns)

    # Decide segments
    if args.segments.strip():
        segs = [s.strip() for s in args.segments.split(",") if s.strip()]
    else:
        segs = [s.strip() for s in (DEFAULT_RSV if platform == "rsv" else DEFAULT_FLU).split(",")]

    metrics = [m.strip() for m in args.keep.split(",") if m.strip()]

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
