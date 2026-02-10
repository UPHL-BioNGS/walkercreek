#!/usr/bin/env python3
import argparse
import pandas as pd

# Fixed, consistent segment order for flu
FLU_SEGMENTS = [
    "A_PB2", "A_PB1", "A_PA", "A_HA", "A_NP", "A_NA", "A_MP", "A_NS",
    "B_PB2", "B_PB1", "B_PA", "B_HA", "B_NP", "B_NA", "B_MP", "B_NS",
]

def normalize_segment(seg: str) -> str:
    """
    Collapse subtype-specific names into stable gene segment names.
    Examples:
      A_HA_H3   -> A_HA
      A_HA_H1   -> A_HA
      A_NA_N2   -> A_NA
      B_HA_*    -> B_HA
      B_NA_*    -> B_NA
    """
    if seg is None:
        return seg
    seg = str(seg).strip()

    # Influenza A
    if seg.startswith("A_HA"):
        return "A_HA"
    if seg.startswith("A_NA"):
        return "A_NA"

    # Influenza B
    if seg.startswith("B_HA"):
        return "B_HA"
    if seg.startswith("B_NA"):
        return "B_NA"

    return seg

def wide_bam(df: pd.DataFrame) -> pd.DataFrame:
    # Expect: Sample, segment_name, number_mapped_reads, mean_depth
    df["segment_norm"] = df["segment_name"].map(normalize_segment)

    # Coerce numerics where possible (keep blanks if not)
    for c in ["number_mapped_reads", "mean_depth"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Build wide tables
    mapped = df.pivot_table(index="Sample", columns="segment_norm", values="number_mapped_reads", aggfunc="first")
    depth  = df.pivot_table(index="Sample", columns="segment_norm", values="mean_depth", aggfunc="first")

    # Force consistent column order, include all segments even if absent
    mapped = mapped.reindex(columns=FLU_SEGMENTS)
    depth  = depth.reindex(columns=FLU_SEGMENTS)

    # Flatten with stable names
    out = pd.DataFrame(index=mapped.index)
    for seg in FLU_SEGMENTS:
        out[f"{seg}_mapped_reads"] = mapped[seg]
        out[f"{seg}_mean_depth"]   = depth[seg]

    out = out.reset_index()
    return out

def wide_coverage(df: pd.DataFrame) -> pd.DataFrame:
    # Expect: Sample, segment_name, reference_length, seq_length, percent_coverage
    df["segment_norm"] = df["segment_name"].map(normalize_segment)

    for c in ["reference_length", "seq_length", "percent_coverage"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    ref = df.pivot_table(index="Sample", columns="segment_norm", values="reference_length", aggfunc="first")
    seqlen = df.pivot_table(index="Sample", columns="segment_norm", values="seq_length", aggfunc="first")
    pcov = df.pivot_table(index="Sample", columns="segment_norm", values="percent_coverage", aggfunc="first")

    ref = ref.reindex(columns=FLU_SEGMENTS)
    seqlen = seqlen.reindex(columns=FLU_SEGMENTS)
    pcov = pcov.reindex(columns=FLU_SEGMENTS)

    out = pd.DataFrame(index=ref.index)
    for seg in FLU_SEGMENTS:
        out[f"{seg}_reference_length"] = ref[seg]
        out[f"{seg}_seq_length"]       = seqlen[seg]
        out[f"{seg}_percent_coverage"] = pcov[seg]

    out = out.reset_index()
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV (normal long format)")
    ap.add_argument("--mode", choices=["bam", "coverage"], required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.inp, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()
    if "Sample" not in df.columns or "segment_name" not in df.columns:
        raise SystemExit(f"Input must contain Sample and segment_name columns. Columns: {list(df.columns)}")

    df["Sample"] = df["Sample"].astype(str).str.strip()

    if args.mode == "bam":
        needed = {"number_mapped_reads", "mean_depth"}
    else:
        needed = {"reference_length", "seq_length", "percent_coverage"}

    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns for mode={args.mode}: {missing}. Columns: {list(df.columns)}")

    out = wide_bam(df) if args.mode == "bam" else wide_coverage(df)

    # Round numeric columns lightly (keep NaN as blank)
    for c in out.columns:
        if c == "Sample":
            continue
        out[c] = pd.to_numeric(out[c], errors="ignore")
    num_cols = out.select_dtypes(include="number").columns
    out[num_cols] = out[num_cols].round(2)

    out.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
