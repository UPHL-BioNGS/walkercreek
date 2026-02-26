#!/usr/bin/env python3

import argparse
import csv
import re
from collections import OrderedDict


def read_tsv(path):
    rows = []
    with open(path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            clean = {str(k).strip(): (v.strip() if isinstance(v, str) else v) for k, v in row.items() if k is not None}
            rows.append(clean)
    return rows


def write_tsv(path, rows, field_order=None):
    if not rows:
        with open(path, "w", newline="", encoding="utf-8") as f:
            f.write("")
        return

    all_fields = []
    seen = set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                seen.add(k)
                all_fields.append(k)

    if field_order:
        ordered = [f for f in field_order if f in seen]
        ordered += [f for f in all_fields if f not in ordered]
    else:
        ordered = all_fields

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=ordered, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def index_by_sample(rows):
    out = OrderedDict()
    for r in rows:
        sample = r.get("Sample", "").strip()
        if not sample:
            continue
        out[sample] = r
    return out


def outer_merge(base_rows, extra_rows):
    base = index_by_sample(base_rows)
    extra = index_by_sample(extra_rows)
    all_samples = list(OrderedDict.fromkeys(list(base.keys()) + list(extra.keys())))
    merged = []
    for sample in all_samples:
        row = OrderedDict()
        row["Sample"] = sample
        if sample in base:
            for k, v in base[sample].items():
                if k != "Sample":
                    row[k] = v
        if sample in extra:
            for k, v in extra[sample].items():
                if k != "Sample" and (k not in row or row[k] == ""):
                    row[k] = v
        merged.append(row)
    return merged


def norm_subtype(raw):
    if raw is None:
        return "NO_SUBTYPE"
    s = str(raw).strip().upper().replace(" ", "").replace("-", "")
    if s == "" or "NOIRMASUBTYPE" in s or "NOABRICATESUBTYPE" in s:
        m_n = re.search(r"(N\d+)", s)
        if m_n and "H" not in s:
            return m_n.group(1)
        return "NO_SUBTYPE"
    if "YAMAGATA" in s or s in {"BYAM", "B/YAM", "YAM"}:
        return "BYAM"
    if "VICTORIA" in s or s in {"BVIC", "B/VIC", "VIC"}:
        return "BVIC"
    m_hn = re.search(r"(H\d+)(N\d+)", s)
    if m_hn:
        return f"{m_hn.group(1)}{m_hn.group(2)}"
    m_h = re.search(r"(H\d+)", s)
    m_n = re.search(r"(N\d+)", s)
    if m_h and m_n:
        return f"{m_h.group(1)}{m_n.group(1)}"
    if m_h:
        return m_h.group(1)
    if m_n:
        return m_n.group(1)
    return "NO_SUBTYPE"


def parse_percent_like(val):
    if val is None:
        return None
    s = str(val).strip()
    if s == "":
        return None
    s = s.replace("%", "").strip()
    try:
        return float(s)
    except ValueError:
        return None


def parse_numeric_prefix(val):
    if val is None:
        return None
    s = str(val).strip()
    if s == "":
        return None
    m = re.search(r"[-+]?\d*\.?\d+", s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None


def is_seasonal_a_subtype(st):
    return st in {"H1N1", "H3N2"}


def is_influenza_b_subtype(st):
    # normalized B lineage in this pipeline
    return st in {"BVIC", "BYAM", "B"}


def is_nonseasonal_a_subtype(st):
    m = re.match(r"^(H\d+)(N\d+)?$", st or "")
    if not m:
        return False
    return m.group(1) not in {"H1", "H3"}


def compute_alerts(row):
    """
    Alerts are intentionally narrow (per Tom's requirements):
        1) Nonseasonal influenza A subtypes (H5/H7/etc) from IRMA
        2) Discordant subtype calls between IRMA and ABRICATE when BOTH are real
        3) Possible coinfection (A + B evidence)
    We do NOT alert on missing subtype/type (RNA degradation noise).
    """
    irma_type = (row.get("IRMA_type") or "").strip()
    irma_sub_raw = row.get("IRMA_subtype", "")
    abr_sub_raw = row.get("abricate_InsaFlu_subtype", "")

    irma_sub = norm_subtype(irma_sub_raw)
    abr_sub = norm_subtype(abr_sub_raw)

    # Evidence for influenza A / B mapping from coverage table
    a_ha_depth = parse_numeric_prefix(row.get("A_HA_mean_depth"))
    a_na_depth = parse_numeric_prefix(row.get("A_NA_mean_depth"))
    a_ha_cov = parse_percent_like(row.get("A_HA_percent_coverage"))
    a_na_cov = parse_percent_like(row.get("A_NA_percent_coverage"))

    b_ha_depth = parse_numeric_prefix(row.get("B_HA_mean_depth"))
    b_na_depth = parse_numeric_prefix(row.get("B_NA_mean_depth"))
    b_ha_cov = parse_percent_like(row.get("B_HA_percent_coverage"))
    b_na_cov = parse_percent_like(row.get("B_NA_percent_coverage"))

    # Kraken %
    flu_a_pct = parse_percent_like(row.get("kraken2 Influenza A percentage"))
    flu_b_pct = parse_percent_like(row.get("kraken2 Influenza B percentage"))

    # Tuneable thresholds for "decent" segment evidence (set at 50% coverage or 50x mean depth to reduce noise)
    def decent_segment(depth, cov):
        return ((cov is not None and cov >= 50.0) or (depth is not None and depth >= 50.0))

    a_evidence = (
        irma_type == "Type_A"
        or (flu_a_pct is not None and flu_a_pct >= 1.0)
        or decent_segment(a_ha_depth, a_ha_cov)
        or decent_segment(a_na_depth, a_na_cov)
    )

    b_evidence = (
        (flu_b_pct is not None and flu_b_pct >= 1.0)
        or decent_segment(b_ha_depth, b_ha_cov)
        or decent_segment(b_na_depth, b_na_cov)
    )

    priority_flag = "NONE"
    discordance_flag = "NONE"
    coinfection_flag = "NONE"
    review_recommended = "NO"

    # 1) Nonseasonal influenza A (from IRMA)
    if irma_type == "Type_A" and is_nonseasonal_a_subtype(irma_sub):
        priority_flag = f"ALERT_NONSEASONAL_A_{irma_sub}"
        review_recommended = "YES"

    # 2) Discordant IRMA vs ABRICATE subtype (only if BOTH are real)
    #    - ignore NO_SUBTYPE on either side
    if irma_sub != "NO_SUBTYPE" and abr_sub != "NO_SUBTYPE" and irma_sub != abr_sub:
        # A vs A discordance (seasonal vs seasonal OR nonseasonal vs seasonal)
        if (is_seasonal_a_subtype(irma_sub) or is_nonseasonal_a_subtype(irma_sub)) and (
            is_seasonal_a_subtype(abr_sub) or is_nonseasonal_a_subtype(abr_sub)
        ):
            discordance_flag = f"ALERT_DISCORDANT_A_SUBTYPE_IRMA_{irma_sub}_ABR_{abr_sub}"
            review_recommended = "YES"
        # A vs B discordance
        elif (is_seasonal_a_subtype(irma_sub) or is_nonseasonal_a_subtype(irma_sub)) and is_influenza_b_subtype(abr_sub):
            discordance_flag = f"ALERT_A_VS_B_SUBTYPE_IRMA_{irma_sub}_ABR_{abr_sub}"
            review_recommended = "YES"
        elif is_influenza_b_subtype(irma_sub) and (is_seasonal_a_subtype(abr_sub) or is_nonseasonal_a_subtype(abr_sub)):
            discordance_flag = f"ALERT_B_VS_A_SUBTYPE_IRMA_{irma_sub}_ABR_{abr_sub}"
            review_recommended = "YES"

    # 3) Possible coinfection: A subtype call + B evidence (or B call + A evidence)
    if (irma_type == "Type_A" and (is_seasonal_a_subtype(irma_sub) or is_nonseasonal_a_subtype(irma_sub))) and b_evidence:
        coinfection_flag = "ALERT_POSSIBLE_COINFECTION_A_PLUS_B"
        review_recommended = "YES"

    if is_influenza_b_subtype(abr_sub) and a_evidence:
        if coinfection_flag == "NONE":
            coinfection_flag = "ALERT_POSSIBLE_COINFECTION_B_PLUS_A"
        review_recommended = "YES"

    # Nextclade status: only meaningful for nonseasonal alerts
    nextclade_qc = (row.get("Nextclade_qc.overallStatus") or "").strip()
    nextclade_seq = (row.get("Nextclade_seqName") or "").strip()
    nextclade_status_flag = "NONE"
    if priority_flag.startswith("ALERT_NONSEASONAL_A_"):
        if nextclade_qc == "" and nextclade_seq == "":
            nextclade_status_flag = "NOT_RUN_NONSEASONAL_SUBTYPE"
        else:
            nextclade_status_flag = "RAN_NONSEASONAL_SUBTYPE"

    # Keep normalized subtype columns (useful for alerts file)
    row["IRMA_subtype_normalized"] = irma_sub
    row["abricate_subtype_normalized"] = abr_sub

    # Alert columns (used for alerts table, later removed from main summary)
    row["priority_subtype_flag"] = priority_flag
    row["subtype_discordance_flag"] = discordance_flag
    row["coinfection_flag"] = coinfection_flag
    row["nextclade_status_flag"] = nextclade_status_flag
    row["review_recommended"] = review_recommended

    return row


def build_alert_rows(rows):
    """Create a compact surveillance alerts table (nonseasonal A + discordance + possible coinfection)."""
    alert_cols = [
        "Sample",
        "IRMA_type",
        "IRMA_subtype",
        "IRMA_subtype_normalized",
        "abricate_InsaFlu_type",
        "abricate_InsaFlu_subtype",
        "abricate_subtype_normalized",

        "priority_subtype_flag",
        "subtype_discordance_flag",
        "coininfection_flag",
        "nextclade_status_flag",
        "review_recommended",

        "Nextclade_qc.overallStatus",
        "Nextclade_seqName",

        "kraken2 Influenza A percentage",
        "kraken2 Influenza B percentage",

        "A_HA_mean_depth",
        "A_HA_percent_coverage",
        "A_NA_mean_depth",
        "A_NA_percent_coverage",

        "B_HA_mean_depth",
        "B_HA_percent_coverage",
        "B_NA_mean_depth",
        "B_NA_percent_coverage",
    ]

    out = []
    for r in rows:
        priority = (r.get("priority_subtype_flag") or "NONE")
        discord = (r.get("subtype_discordance_flag") or "NONE")
        coinfect = (r.get("coininfection_flag") or "NONE")
        review = (r.get("review_recommended") or "NO")

        if priority != "NONE" or discord != "NONE" or coinfect != "NONE" or review == "YES":
            row = OrderedDict()
            for c in alert_cols:
                row[c] = r.get(c, "")
            out.append(row)

    return out, alert_cols


def strip_alert_columns_for_summary(rows):
    """Remove helper/alert columns from the main summary report."""
    cols_to_remove = {
        "IRMA_subtype_normalized",
        "abricate_subtype_normalized",
        "priority_subtype_flag",
        "subtype_discordance_flag",
        "coininfection_flag",
        "nextclade_status_flag",
        "review_recommended",
    }

    cleaned = []
    for r in rows:
        row = OrderedDict((k, v) for k, v in r.items() if k not in cols_to_remove)
        cleaned.append(row)
    return cleaned


def main():
    parser = argparse.ArgumentParser(
        description="Merge Walkercreek report TSVs and write a compact surveillance alerts table."
    )

    parser.add_argument("--qc", required=True, help="QC reportsheet TSV")
    parser.add_argument("--typing", required=True, help="IRMA+ABRICATE typing report TSV")
    parser.add_argument("--irma-qc", required=True, dest="irma_qc", help="IRMA consensus QC TSV")
    parser.add_argument("--nextclade", required=True, help="Nextclade report TSV")
    parser.add_argument("--coverage", required=True, help="Merged BAM+coverage TSV")
    parser.add_argument("--kraken2", required=False, help="Kraken2 reportsheet TSV")
    parser.add_argument("--out", default="summary_report.tsv", help="Output summary TSV")
    parser.add_argument("--alerts-out", default="summary_alerts.tsv", help="Output surveillance alerts TSV")

    # Backward-compatible positional args support
    parser.add_argument("legacy", nargs="*", help=argparse.SUPPRESS)

    args = parser.parse_args()

    # Support old positional calling style:
    # with kraken:  qc typing irma_qc nextclade kraken2 coverage
    # no kraken:    qc typing nextclade irma_qc coverage   (older order)
    if args.legacy:
        n = len(args.legacy)
        if n == 6:
            args.qc = args.legacy[0]
            args.typing = args.legacy[1]
            args.irma_qc = args.legacy[2]
            args.nextclade = args.legacy[3]
            args.kraken2 = args.legacy[4]
            args.coverage = args.legacy[5]
        elif n == 5:
            args.qc = args.legacy[0]
            args.typing = args.legacy[1]
            args.nextclade = args.legacy[2]
            args.irma_qc = args.legacy[3]
            args.coverage = args.legacy[4]
        else:
            raise SystemExit(f"Unsupported positional argument count: {n}. Expected 5 or 6.")

    qc_rows = read_tsv(args.qc)
    typing_rows = read_tsv(args.typing)
    irma_qc_rows = read_tsv(args.irma_qc)
    nextclade_rows = read_tsv(args.nextclade)
    coverage_rows = read_tsv(args.coverage)
    kraken_rows = read_tsv(args.kraken2) if args.kraken2 else []

    merged = qc_rows
    for extra in [typing_rows, irma_qc_rows, nextclade_rows, coverage_rows]:
        merged = outer_merge(merged, extra)

    if kraken_rows:
        merged = outer_merge(merged, kraken_rows)

    merged = [compute_alerts(r) for r in merged]

    # Write compact alerts table
    alert_rows, alert_order = build_alert_rows(merged)
    write_tsv(args.alerts_out, alert_rows, field_order=alert_order)

    # Remove alert/helper columns from main summary report
    summary_rows = strip_alert_columns_for_summary(merged)

    preferred_order = [
        "Sample",
        "reads_before_trimming",
        "GC_before_trimming",
        "average_Q_score_before_trimming",
        "reads_after_trimming",
        "paired_reads_after_trimming",
        "unpaired_reads_after_trimming",
        "GC_after_trimming",
        "average_Q_score_after_trimming",

        "IRMA_type",
        "IRMA_subtype",
        "abricate_InsaFlu_type",
        "abricate_InsaFlu_subtype",

        "IRMA_consensus_ACTG_count",
        "IRMA_consensus_degenerate_count",
        "IRMA_consensus_N_count",
        "IRMA_consensus_total_count",
        "IRMA_consensus_segment_count",
        "IRMA_consensus_N50",
        "IRMA_consensus_GC_content",

        "clade",
        "legacy_clade",
        "short_clade",
        "subclade",
        "Nextclade_qc.overallStatus",
        "Nextclade_totalSubstitutions",
        "Nextclade_coverage",
        "Nextclade_seqName",

        "kraken2 Homo sapiens percentage",
        "kraken2 Influenza A percentage",
        "kraken2 Influenza B percentage",
        "kraken2 unclassified percentage",
    ]

    write_tsv(args.out, summary_rows, field_order=preferred_order)


if __name__ == "__main__":
    main()
