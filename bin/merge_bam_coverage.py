#!/usr/bin/env python3
import argparse
import sys
import pandas as pd

# segment order for FLU
FLU_SEGMENTS = [
    "A_PB2",
    "A_PB1",
    "A_PA",
    "A_HA",
    "A_NP",
    "A_NA",
    "A_MP",
    "A_NS",
    "B_PB2",
    "B_PB1",
    "B_PA",
    "B_HA",
    "B_NP",
    "B_NA",
    "B_MP",
    "B_NS",
]

# segment order for RSV
RSV_SEGMENTS = ["RSV_A", "RSV_AD", "RSV_B", "RSV_BD"]

BAM_LONG_NEEDED = {"Sample", "segment_name", "number_mapped_reads", "mean_depth"}
COV_LONG_NEEDED = {"Sample", "segment_name", "reference_length", "seq_length", "percent_coverage"}


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
    return BAM_LONG_NEEDED.issubset(set(df.columns))


def is_long_cov(df: pd.DataFrame) -> bool:
    return COV_LONG_NEEDED.issubset(set(df.columns))


def detect_platform_from_segments(segments: list[str]) -> str:
    """
    Return 'rsv' or 'flu' based on segment name patterns.
    """
    segs = [str(s).strip() for s in segments if s is not None]
    # RSV segments look like RSV_A / RSV_AD / RSV_B / RSV_BD
    if any(s.startswith("RSV_") for s in segs):
        return "rsv"
    # FLU segments look like A_* / B_*
    if any(s.startswith("A_") or s.startswith("B_") for s in segs):
        return "flu"
    # fallback
    return "flu"


def normalize_segment(seg: str, platform: str) -> str:
    if seg is None:
        return seg
    seg = str(seg).strip()

    # For RSV, keep exactly as observed (RSV_A/RSV_AD/RSV_B/RSV_BD)
    if platform == "rsv":
        return seg

    # FLU normalization
    if seg.startswith("A_HA"):
        return "A_HA"
    if seg.startswith("A_NA"):
        return "A_NA"
    if seg.startswith("B_HA"):
        return "B_HA"
    if seg.startswith("B_NA"):
        return "B_NA"
    return seg


def long_bam_to_wide(df: pd.DataFrame, segment_order: list[str], platform: str) -> pd.DataFrame:
    df = df.copy()
    df["segment_norm"] = df["segment_name"].map(lambda s: normalize_segment(s, platform))

    for c in ["number_mapped_reads", "mean_depth"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    mapped = df.pivot_table(index="Sample", columns="segment_norm", values="number_mapped_reads", aggfunc="first")
    depth = df.pivot_table(index="Sample", columns="segment_norm", values="mean_depth", aggfunc="first")

    mapped = mapped.reindex(columns=segment_order)
    depth = depth.reindex(columns=segment_order)

    out = pd.DataFrame(index=mapped.index)
    for seg in segment_order:
        out[f"{seg}_mapped_reads"] = mapped[seg]
        out[f"{seg}_mean_depth"] = depth[seg]
    out = out.reset_index()
    return out


def long_cov_to_wide(df: pd.DataFrame, segment_order: list[str], platform: str) -> pd.DataFrame:
    df = df.copy()
    df["segment_norm"] = df["segment_name"].map(lambda s: normalize_segment(s, platform))

    for c in ["reference_length", "seq_length", "percent_coverage"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    ref = df.pivot_table(index="Sample", columns="segment_norm", values="reference_length", aggfunc="first")
    seqlen = df.pivot_table(index="Sample", columns="segment_norm", values="seq_length", aggfunc="first")
    pcov = df.pivot_table(index="Sample", columns="segment_norm", values="percent_coverage", aggfunc="first")

    ref = ref.reindex(columns=segment_order)
    seqlen = seqlen.reindex(columns=segment_order)
    pcov = pcov.reindex(columns=segment_order)

    out = pd.DataFrame(index=ref.index)
    for seg in segment_order:
        out[f"{seg}_percent_coverage"] = pcov[seg]
        out[f"{seg}_reference_length"] = ref[seg]
        out[f"{seg}_seq_length"] = seqlen[seg]
    out = out.reset_index()
    return out


def force_all_columns(df: pd.DataFrame, segment_order: list[str]) -> pd.DataFrame:
    """
    Ensure df has every expected column in canonical order.
    Missing columns are added as NA.
    """
    expected = ["Sample"]
    for seg in segment_order:
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
    ap.add_argument(
        "--platform",
        default="auto",
        choices=["auto", "flu", "rsv"],
        help="Force platform schema. Default auto-detect from segment_name values.",
    )
    args = ap.parse_args()

    bam = read_tsv(args.bam)
    cov = read_tsv(args.cov)

    platform = args.platform

    segs = []

    if platform == "auto":
        if "segment_name" in bam.columns:
            segs.extend(bam["segment_name"].dropna().tolist())
        if "segment_name" in cov.columns:
            segs.extend(cov["segment_name"].dropna().tolist())

        if segs:
            platform = detect_platform_from_segments(segs)
        else:
            # Wide-table fallback: detect by column prefixes
            cols = list(bam.columns) + list(cov.columns)
            if any(str(c).startswith("RSV_") for c in cols):
                platform = "rsv"
            elif any(str(c).startswith("A_") or str(c).startswith("B_") for c in cols):
                platform = "flu"
            else:
                platform = "flu"

    segment_order = RSV_SEGMENTS if platform == "rsv" else FLU_SEGMENTS

    # Convert long->wide if needed
    bam_w = long_bam_to_wide(bam, segment_order, platform) if is_long_bam(bam) else bam
    cov_w = long_cov_to_wide(cov, segment_order, platform) if is_long_cov(cov) else cov

    # Merge on Sample
    merged = pd.merge(bam_w, cov_w, on="Sample", how="outer")

    # Normalize ordering + ensure all columns exist
    merged = force_all_columns(merged, segment_order)
    merged = coerce_numeric_round(merged)

    merged.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
