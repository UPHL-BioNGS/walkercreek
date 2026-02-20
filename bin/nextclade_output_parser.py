#!/usr/bin/env python3

import argparse
import pandas as pd


def first_present(row_df: pd.DataFrame, candidates, default=""):
    """
    Return the first column found in row_df from `candidates` as a scalar.
    If none exist, return the default.
    """
    for c in candidates:
        if c in row_df.columns:
            v = row_df[c].iloc[0]
            # Normalize NaN -> ""
            if pd.isna(v):
                return default
            return v
    return default


def as_float_or_blank(v):
    try:
        if v == "" or v is None or (isinstance(v, float) and pd.isna(v)):
            return ""
        return float(v)
    except Exception:
        return ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--id", required=True, help="Path to Nextclade TSV output (e.g. sample.tsv)")
    args = parser.parse_args()

    tsv_path = args.id
    sample_name = tsv_path.rsplit(".", 1)[0]  # remove .tsv

    # Read TSV
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as e:
        raise SystemExit(f"ERROR: Could not read Nextclade TSV '{tsv_path}': {e}")

    if df.empty:
        raise SystemExit(f"ERROR: Nextclade TSV '{tsv_path}' is empty")

    clade_val = ""
    if "clade" in df.columns and not df["clade"].empty and pd.notna(df["clade"].iloc[0]):
        clade_val = str(df["clade"].iloc[0])

    # Always create the file (even if blank)
    with open("NEXTCLADE_CLADE.tsv", "w") as fh:
        fh.write(f"0\t{clade_val}\n")

    # Nextclade is typically 1 row per sequence; take row 0
    row = df.iloc[[0]].copy()

    # Convert coverage fraction -> percent if present
    if "coverage" in row.columns:
        try:
            row["coverage"] = (row["coverage"].astype(float) * 100).round(2)
        except Exception:
            pass  # leave as-is if unexpected type

    # Build summary
    out = {
        "Sample": sample_name,
        # Core clade fields
        "clade": first_present(row, ["clade"], default=""),
        "legacy_clade": first_present(row, ["legacy-clade", "legacy_clade"], default=""),
        "short_clade": first_present(row, ["short-clade", "short_clade"], default=""),
        "subclade": first_present(row, ["subclade"], default=""),
        # Core QC
        "Nextclade_qc.overallStatus": first_present(row, ["qc.overallStatus"], default=""),
        # Minimal counts
        "Nextclade_totalSubstitutions": first_present(row, ["totalSubstitutions"], default=""),
        "Nextclade_coverage": first_present(row, ["coverage"], default=""),
        # High-value debug
        "Nextclade_seqName": first_present(row, ["seqName"], default=""),
    }

    # Coerce numeric fields that often come out as floats/ints to keep TSV clean (leave blanks if missing)
    for k in ["Nextclade_totalSubstitutions", "Nextclade_coverage"]:
        out[k] = out[k] if out[k] == "" else out[k]

    out_df = pd.DataFrame([out])

    # Optionally round numeric columns if they parse as numbers
    for col in ["Nextclade_coverage"]:
        try:
            out_df[col] = pd.to_numeric(out_df[col], errors="ignore")
            if pd.api.types.is_numeric_dtype(out_df[col]):
                out_df[col] = out_df[col].round(2)
        except Exception:
            pass

    out_df.to_csv(f"{sample_name}.nextclade_report.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
