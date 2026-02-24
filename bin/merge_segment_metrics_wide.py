#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

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

RSV_SEGMENTS = ["RSV_A", "RSV_AD", "RSV_B", "RSV_BD"]


def detect_platform_from_segments(segs) -> str:
    segs = [str(s).strip() for s in segs if s is not None]
    if any(s.startswith("RSV_") for s in segs):
        return "rsv"
    if any(s.startswith("A_") or s.startswith("B_") for s in segs):
        return "flu"
    return "flu"


def normalize_segment(seg: str, platform: str) -> str:
    """
    FLU: collapse subtype-specific names into stable gene segment names.
    RSV: keep segment_name as-is (RSV_A/RSV_AD/RSV_B/RSV_BD).
    """
    if seg is None:
        return seg
    seg = str(seg).strip()

    if platform == "rsv":
        return seg

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


def wide_bam(df: pd.DataFrame, segment_order, platform: str) -> pd.DataFrame:
    # Expect: Sample, segment_name, number_mapped_reads, mean_depth
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

    return out.reset_index()


def wide_coverage(df: pd.DataFrame, segment_order, platform: str) -> pd.DataFrame:
    # Expect: Sample, segment_name, reference_length, seq_length, percent_coverage
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
        out[f"{seg}_reference_length"] = ref[seg]
        out[f"{seg}_seq_length"] = seqlen[seg]
        out[f"{seg}_percent_coverage"] = pcov[seg]

    return out.reset_index()


def round_numeric(out: pd.DataFrame) -> pd.DataFrame:
    out = out.copy()
    for c in out.columns:
        if c == "Sample":
            continue
        out[c] = pd.to_numeric(out[c], errors="coerce")
    num_cols = out.select_dtypes(include="number").columns
    out[num_cols] = out[num_cols].round(2)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV (long format)")
    ap.add_argument("--mode", choices=["bam", "coverage"], required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--platform",
        choices=["auto", "flu", "rsv"],
        default="auto",
        help="Platform schema. Default auto-detect from segment_name.",
    )
    args = ap.parse_args()

    df = pd.read_csv(args.inp, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()

    if "Sample" not in df.columns or "segment_name" not in df.columns:
        raise SystemExit(f"Input must contain Sample and segment_name. Columns: {list(df.columns)}")

    df["Sample"] = df["Sample"].astype(str).str.strip()

    platform = args.platform
    if platform == "auto":
        platform = detect_platform_from_segments(df["segment_name"].dropna().tolist())

    segment_order = RSV_SEGMENTS if platform == "rsv" else FLU_SEGMENTS

    if args.mode == "bam":
        needed = {"number_mapped_reads", "mean_depth"}
        missing = [c for c in needed if c not in df.columns]
        if missing:
            raise SystemExit(f"Missing required columns for mode=bam: {missing}. Columns: {list(df.columns)}")
        out = wide_bam(df, segment_order, platform)

    else:
        needed = {"reference_length", "seq_length", "percent_coverage"}
        missing = [c for c in needed if c not in df.columns]
        if missing:
            raise SystemExit(f"Missing required columns for mode=coverage: {missing}. Columns: {list(df.columns)}")
        out = wide_coverage(df, segment_order, platform)

    out = round_numeric(out)
    out.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
