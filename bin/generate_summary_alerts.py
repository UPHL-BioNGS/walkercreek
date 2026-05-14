#!/usr/bin/env python3

"""
Generate summary_alerts.tsv for walkercreek summary processes.

This script is intentionally defensive:
- always writes a header, even if no alerts are found
- scans TSV files in the current work directory unless specific files are provided
- groups evidence by sample ID
- detects likely influenza/RSV mixed detections and coinfections
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

HEADER = [
    "sample",
    "alert_type",
    "severity",
    "message",
    "evidence",
    "source_files",
]


SAMPLE_COL_CANDIDATES = [
    "sample",
    "sample_id",
    "sampleid",
    "id",
    "sample_name",
    "sample_name_alias",
    "name",
    "run_accession",
    "sra",
    "accession",
]


# Specific detections are more useful than generic "Flu A"/"Flu B"
SPECIFIC_TOKEN_PATTERNS = {
    "H1N1": re.compile(r"\bH1N1\b|A/H1|H1N1pdm|pdm09", re.I),
    "H3N2": re.compile(r"\bH3N2\b|A/H3", re.I),
    "H5NX": re.compile(r"\bH5N[X0-9]+\b|\bH5\b|H5NX", re.I),
    "H7NX": re.compile(r"\bH7N[X0-9]+\b|\bH7\b", re.I),
    "B_VIC": re.compile(r"\bB[_\-\s]?VIC\b|Victoria", re.I),
    "B_YAM": re.compile(r"\bB[_\-\s]?YAM\b|Yamagata", re.I),
    "RSV_A": re.compile(r"\bRSV[_\-\s]?A\b|\bRSV-A\b", re.I),
    "RSV_B": re.compile(r"\bRSV[_\-\s]?B\b|\bRSV-B\b", re.I),
}

GENERIC_TOKEN_PATTERNS = {
    "FLU_A": re.compile(r"\bInfluenza\s*A\b|\bFlu\s*A\b|\bType\s*A\b", re.I),
    "FLU_B": re.compile(r"\bInfluenza\s*B\b|\bFlu\s*B\b|\bType\s*B\b", re.I),
    "RSV": re.compile(r"\bRSV\b|respiratory\s+syncytial", re.I),
}


NEGATIVE_PATTERNS = [
    re.compile(r"\bnot\s+detected\b", re.I),
    re.compile(r"\bnegative\b", re.I),
    re.compile(r"\bnone\b", re.I),
    re.compile(r"\bno\s+hit\b", re.I),
    re.compile(r"\bno\s+subtype\b", re.I),
    re.compile(r"\bunknown\b", re.I),
    re.compile(r"\bNA\b"),
    re.compile(r"^\s*$"),
]


def is_negative(value: str) -> bool:
    text = str(value or "").strip()
    return any(p.search(text) for p in NEGATIVE_PATTERNS)


def normalize_header(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def sniff_sample_column(fieldnames: List[str]) -> str | None:
    normalized = {normalize_header(c): c for c in fieldnames}
    for candidate in SAMPLE_COL_CANDIDATES:
        if candidate in normalized:
            return normalized[candidate]

    # fallback: first column with "sample" in name
    for c in fieldnames:
        if "sample" in normalize_header(c):
            return c

    return fieldnames[0] if fieldnames else None


def read_tsv(path: Path) -> Iterable[Tuple[str, Dict[str, str]]]:
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                return

            sample_col = sniff_sample_column(reader.fieldnames)
            if not sample_col:
                return

            for row in reader:
                sample = str(row.get(sample_col, "")).strip()
                if not sample:
                    continue
                yield sample, row

    except UnicodeDecodeError:
        return
    except Exception:
        return


def tokens_from_row(row: Dict[str, str]) -> Set[str]:
    tokens: Set[str] = set()

    joined = " ".join(str(v or "") for v in row.values())

    # Skip rows that are entirely negative/empty.
    if is_negative(joined):
        return tokens

    for token, pattern in SPECIFIC_TOKEN_PATTERNS.items():
        if pattern.search(joined):
            tokens.add(token)

    # Add generic tokens only if useful.
    for token, pattern in GENERIC_TOKEN_PATTERNS.items():
        if pattern.search(joined):
            tokens.add(token)

    return tokens


def simplify_tokens(tokens: Set[str]) -> Set[str]:
    simplified = set(tokens)

    # If specific Flu A subtypes exist, generic FLU_A is not needed.
    if {"H1N1", "H3N2", "H5NX", "H7NX"} & simplified:
        simplified.discard("FLU_A")

    # If specific Flu B lineage exists, generic FLU_B is not needed.
    if {"B_VIC", "B_YAM"} & simplified:
        simplified.discard("FLU_B")

    # If specific RSV subtypes exist, generic RSV is not needed.
    if {"RSV_A", "RSV_B"} & simplified:
        simplified.discard("RSV")

    return simplified


def classify_alert(tokens: Set[str]) -> Tuple[str | None, str | None, str | None]:
    flu_a = tokens & {"H1N1", "H3N2", "H5NX", "H7NX", "FLU_A"}
    flu_b = tokens & {"B_VIC", "B_YAM", "FLU_B"}
    rsv = tokens & {"RSV_A", "RSV_B", "RSV"}

    if len(flu_a) >= 2:
        return (
            "influenza_a_mixed_subtype",
            "warning",
            f"Multiple influenza A subtype signals detected: {', '.join(sorted(flu_a))}",
        )

    if len(flu_b) >= 2:
        return (
            "influenza_b_mixed_lineage",
            "warning",
            f"Multiple influenza B lineage signals detected: {', '.join(sorted(flu_b))}",
        )

    if flu_a and flu_b:
        return (
            "influenza_a_b_coinfection",
            "warning",
            f"Influenza A and influenza B signals detected: {', '.join(sorted(flu_a | flu_b))}",
        )

    if len(rsv) >= 2:
        return (
            "rsv_mixed_subtype",
            "warning",
            f"Multiple RSV subtype signals detected: {', '.join(sorted(rsv))}",
        )

    # Broader mixed virus signal. Useful but lower confidence.
    virus_groups = sum(bool(x) for x in [flu_a, flu_b, rsv])
    if virus_groups >= 2:
        return (
            "mixed_virus_signal",
            "info",
            f"Multiple respiratory virus groups detected: {', '.join(sorted(tokens))}",
        )

    return None, None, None


def collect_files(paths: List[str]) -> List[Path]:
    if paths:
        files = [Path(p) for p in paths]
    else:
        files = sorted(Path(".").glob("*.tsv"))

    ignored_names = {
        "summary_alerts.tsv",
        "versions.yml",
    }

    return [p for p in files if p.is_file() and p.name not in ignored_names and p.suffix.lower() == ".tsv"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="summary_alerts.tsv")
    parser.add_argument("files", nargs="*")
    args = parser.parse_args()

    files = collect_files(args.files)

    evidence_by_sample: Dict[str, Set[str]] = defaultdict(set)
    source_by_sample: Dict[str, Set[str]] = defaultdict(set)

    for path in files:
        for sample, row in read_tsv(path):
            tokens = tokens_from_row(row)
            if not tokens:
                continue

            tokens = simplify_tokens(tokens)
            if not tokens:
                continue

            evidence_by_sample[sample].update(tokens)
            source_by_sample[sample].add(path.name)

    rows: List[List[str]] = []

    for sample in sorted(evidence_by_sample):
        tokens = simplify_tokens(evidence_by_sample[sample])
        alert_type, severity, message = classify_alert(tokens)

        if not alert_type:
            continue

        rows.append(
            [
                sample,
                alert_type,
                severity,
                message,
                ",".join(sorted(tokens)),
                ",".join(sorted(source_by_sample[sample])),
            ]
        )

    out = Path(args.out)
    with out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
