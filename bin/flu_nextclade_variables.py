#!/usr/bin/env python3

import argparse
import csv
import os
import re
from pathlib import Path
from typing import Tuple, List


SEASONAL_DATASETS = {
    "H1N1": "flu_h1n1pdm_ha",
    "H3N2": "flu_h3n2_ha",
    "BVIC": "flu_vic_ha",
    "BYAM": "flu_yam_ha",
}

# Add non-seasonal dataset routing names 
NONSEASONAL_DATASETS = {
    "H5": "flu_h5_ha",
    "H7": "flu_h7_ha",
    "H9": "flu_h9_ha",
    "H10": "flu_h10_ha",
}


def read_first_line(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.readline().strip()


def normalize_subtype(raw: str) -> str:
    """
    Normalize strings like:
        H3N2, h3n2, Type_A H5N1, No IRMA subtype, NoIRMAsubtypeN2, etc.
    Returns normalized values like:
        H3N2, H5N1, H5, N1, BVIC, BYAM, NO_SUBTYPE
    """
    if not raw:
        return "NO_SUBTYPE"

    s = raw.strip().upper().replace(" ", "")
    s = s.replace("-", "")

    null_markers = {
        "", "NOIRMASUBTYPE", "NOIRMASTUBTYPE", "NO_SUBTYPE",
        "NOABRICATESUBTYPE", "NOIRMATYPE",
        "NOABRICATETYPE", "NOIRMA", "NONE", "NA"
    }
    if s in null_markers or "NOIRMASUBTYPE" in s or "NOABRICATESUBTYPE" in s:
        m_n = re.search(r"(N\d+)", s)
        if m_n and "H" not in s:
            return m_n.group(1)
        return "NO_SUBTYPE"

    # Influenza B lineage normalization
    if "YAMAGATA" in s or s in {"BYAM", "B/YAM", "YAM"}:
        return "BYAM"
    if "VICTORIA" in s or s in {"BVIC", "B/VIC", "VIC"}:
        return "BVIC"

    # Extract HxNy / Hx / Ny
    m_hn = re.search(r"(H\d+)(N\d+)", s)
    if m_hn:
        return "{}{}".format(m_hn.group(1), m_hn.group(2))

    m_h = re.search(r"(H\d+)", s)
    m_n = re.search(r"(N\d+)", s)

    if m_h and m_n:
        return "{}{}".format(m_h.group(1), m_n.group(1))
    if m_h:
        return m_h.group(1)
    if m_n:
        return m_n.group(1)

    return "NO_SUBTYPE"


def choose_subtype(irma_subtype: str, abricate_subtype: str) -> Tuple[str, str]:
    """
    Prefer IRMA if usable, fallback to ABRICATE, otherwise NO_SUBTYPE.
    Returns (chosen_subtype, source_used)
    """
    irma_norm = normalize_subtype(irma_subtype)
    abr_norm = normalize_subtype(abricate_subtype)

    if irma_norm != "NO_SUBTYPE":
        return irma_norm, "IRMA"
    if abr_norm != "NO_SUBTYPE":
        return abr_norm, "ABRICATE"
    return "NO_SUBTYPE", "NONE"


def classify_route(chosen_subtype: str) -> Tuple[str, str]:
    """
    Returns (route_type, route_value)
    route_type:
        - seasonal_nextclade
        - nonseasonal_nextclade
        - unresolved_subtype
        - no_subtype
    """
    # Seasonal A
    if chosen_subtype == "H1N1":
        return "seasonal_nextclade", SEASONAL_DATASETS["H1N1"]
    if chosen_subtype == "H3N2":
        return "seasonal_nextclade", SEASONAL_DATASETS["H3N2"]

    # Seasonal B
    if chosen_subtype == "BVIC":
        return "seasonal_nextclade", SEASONAL_DATASETS["BVIC"]
    if chosen_subtype == "BYAM":
        return "seasonal_nextclade", SEASONAL_DATASETS["BYAM"]

    # Non-seasonal influenza A (H5, H7, etc.)
    m_h = re.match(r"^(H\d+)(N\d+)?$", chosen_subtype)
    if m_h:
        ha = m_h.group(1)
        if ha in NONSEASONAL_DATASETS:
            return "nonseasonal_nextclade", NONSEASONAL_DATASETS[ha]
        if ha not in {"H1", "H3"}:
            # Unknown non-seasonal HA: still flag, no dataset route
            return "unresolved_subtype", "nonseasonal_influenza_A_no_dataset:{}".format(chosen_subtype)

    # NA-only or odd subtype remnants
    if chosen_subtype.startswith("N"):
        return "unresolved_subtype", "NA_only:{}".format(chosen_subtype)

    if chosen_subtype == "NO_SUBTYPE":
        return "no_subtype", "no_subtype"

    return "unresolved_subtype", "unrecognized:{}".format(chosen_subtype)


def write_tsv(path: str, fieldnames: List[str], row: dict) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Route flu samples to Nextclade datasets and flag non-seasonal subtypes.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--irma_subtype", required=True)
    parser.add_argument("--abricate_subtype", required=True)
    args = parser.parse_args()

    sample = args.sample
    irma_raw = read_first_line(args.irma_subtype) if os.path.exists(args.irma_subtype) else ""
    abr_raw = read_first_line(args.abricate_subtype) if os.path.exists(args.abricate_subtype) else ""

    chosen_subtype, source_used = choose_subtype(irma_raw, abr_raw)
    route_type, route_value = classify_route(chosen_subtype)

    # Always write route decision file
    route_tsv = "{}.nextclade_route.tsv".format(sample)
    write_tsv(
        route_tsv,
        [
            "Sample",
            "IRMA_subtype_raw",
            "abricate_subtype_raw",
            "chosen_subtype",
            "chosen_subtype_source",
            "nextclade_route_type",
            "nextclade_route_value",
        ],
        {
            "Sample": sample,
            "IRMA_subtype_raw": irma_raw or "NA",
            "abricate_subtype_raw": abr_raw or "NA",
            "chosen_subtype": chosen_subtype,
            "chosen_subtype_source": source_used,
            "nextclade_route_type": route_type,
            "nextclade_route_value": route_value,
        },
    )

    # Route to dataset trigger directories (seasonal OR nonseasonal)
    if route_type in {"seasonal_nextclade", "nonseasonal_nextclade"}:
        Path(route_value).mkdir(parents=True, exist_ok=True)

    # Always write a flag for nonseasonal/unresolved/no-subtype so they stand out in summary
    if route_type in {"nonseasonal_nextclade", "unresolved_subtype", "no_subtype"}:
        if route_type == "nonseasonal_nextclade":
            flag_reason = "nonseasonal_influenza_A:{}".format(chosen_subtype)
        else:
            flag_reason = route_value

        flag_tsv = "{}.nonseasonal_flag.tsv".format(sample)
        write_tsv(
            flag_tsv,
            [
                "Sample",
                "flag_type",
                "flag_reason",
                "chosen_subtype",
                "chosen_subtype_source",
                "IRMA_subtype_raw",
                "abricate_subtype_raw",
                "nextclade_route_type",
                "nextclade_route_value",
            ],
            {
                "Sample": sample,
                "flag_type": route_type,
                "flag_reason": flag_reason,
                "chosen_subtype": chosen_subtype,
                "chosen_subtype_source": source_used,
                "IRMA_subtype_raw": irma_raw or "NA",
                "abricate_subtype_raw": abr_raw or "NA",
                "nextclade_route_type": route_type,
                "nextclade_route_value": route_value,
            },
        )

if __name__ == "__main__":
    main()

