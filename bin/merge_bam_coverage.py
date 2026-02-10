#!/usr/bin/env python3
import argparse
import sys
import pandas as pd

# Fixed, consistent segment order for flu
FLU_SEGMENTS = [
    "A_PB2", "A_PB1", "A_PA", "A_HA", "A_NP", "A_NA", "A_MP", "A_NS",
    "B_PB2", "B_PB1", "B_PA", "B_HA", "B_NP", "B_NA", "B_MP", "B_NS",
]

def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)

def read_tsv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()
    if "Sample" not in df.columns:
        die(f"TSV missing Sample column: {path}. Have={list(df.columns)}")
    df["Sample"] = df["Sample"].astype(str).str.strip()
    return df

def is_long_bam(df: pd.DataFrame) -> bool:
    needed = {"Sample", "segment_name", "number_mapped_reads", "mean_depth"}
    return needed.issubset(set(df.columns))

def is_long_cov(df: pd.DataFrame) -> bool:
    needed = {"Sample", "segment_name", "reference_length", "seq_length", "percent_coverage"}
    return needed.issubset(set(df.columns))

def normalize_segment(seg: str) -> str:
    if seg is None:
        return seg
    seg = str(seg).strip()
    if seg.startswith("A_HA"):
        return "A_HA"
    if seg.startswith("A_NA"):
        return "A_NA"
    if seg.startswith("B_HA"):
        return "B_HA"
    if seg.startswith("B_NA"):
        return "B_NA"
    return seg

def long_bam_to_wide(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["segment_norm"] = df["segment_name"].map(normalize_segment)

    for c in ["number_mapped_reads", "mean_depth"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    mapped = df.pivot_table(index="Sample", columns="segment_norm", values="number_mapped_reads", aggfunc="first")
    depth  = df.pivot_table(index="Sample", columns="segment_norm", values="mean_depth", aggfunc="first")

    mapped = mapped.reindex(columns=FLU_SEGMENTS)
    depth  = depth.reindex(columns=FLU_SEGMENTS)

    out = pd.DataFrame(index=mapped.index)
    for seg in FLU_SEGMENTS:
        out[f"{seg}_mapped_reads"] = mapped[seg]
        out[f"{seg}_mean_depth"]   = depth[seg]
    out = out.reset_index()
    return out

def long_cov_to_wide(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["segment_norm"] = df["segment_name"].map(normalize_segment)

    for c in ["reference_length", "seq_length", "percent_coverage"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    ref    = df.pivot_table(index="Sample", columns="segment_norm", values="reference_length", aggfunc="first")
    seqlen = df.pivot_table(index="Sample", columns="segment_norm", values="seq_length", aggfunc="first")
    pcov   = df.pivot_table(index="Sample", columns="segment_norm", values="percent_coverage", aggfunc="first")

    ref    = ref.reindex(columns=FLU_SEGMENTS)
    seqlen = seqlen.reindex(columns=FLU_SEGMENTS)
    pcov   = pcov.reindex(columns=FLU_SEGMENTS)

    out = pd.DataFrame(index=ref.index)
    for seg in FLU_SEGMENTS:
        out[f"{seg}_percent_coverage"] = pcov[seg]
        out[f"{seg}_reference_length"] = ref[seg]
        out[f"{seg}_seq_length"]       = seqlen[seg]
    out = out.reset_index()
    return out

def force_all_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure df has every expected column in the canonical order.
    Missing columns are added as NA.
    """
    expected = ["Sample"]
    for seg in FLU_SEGMENTS:
        expected += [
            f"{seg}_mapped_reads",
            f"{seg}_mean_depth",
            f"{seg}_percent_coverage",
            f"{seg}_reference_length",
            f"{seg}_seq_length",
        ]

    df = df.copy()
    for col in expected:
        if col not in df.columns:
            df[col] = pd.NA
    df = df[expected]
    return df

def coerce_numeric_round(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for c in df.columns:
        if c == "Sample":
            continue
        df[c] = pd.to_numeric(df[c], errors="coerce")
    num_cols = df.select_dtypes(include="number").columns
    df[num_cols] = df[num_cols].round(2)
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True, help="BAM metrics TSV (wide or long)")
    ap.add_argument("--cov", required=True, help="Coverage metrics TSV (wide or long)")
    ap.add_argument("--out", required=True, help="Output TSV")
    args = ap.parse_args()

    bam = read_tsv(args.bam)
    cov = read_tsv(args.cov)

    # Convert long->wide if needed
    if is_long_bam(bam):
        bam_w = long_bam_to_wide(bam)
    else:
        bam_w = bam

    if is_long_cov(cov):
        cov_w = long_cov_to_wide(cov)
    else:
        cov_w = cov

    # Merge
    merged = pd.merge(bam_w, cov_w, on="Sample", how="outer")

    # Normalize ordering + ensure all segment columns exist
    merged = force_all_columns(merged)
    merged = coerce_numeric_round(merged)

    merged.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
