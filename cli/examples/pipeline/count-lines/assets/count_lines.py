#!/usr/bin/env python3
"""Annotate a samplesheet with per-file line counts."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to the input samplesheet CSV")
    parser.add_argument("--output", required=True, help="Path to write the annotated CSV")
    parser.add_argument(
        "--assets-dir",
        required=True,
        help="Directory containing files referenced by the samplesheet",
    )
    return parser.parse_args()


def detect_file_column(fieldnames: Iterable[str]) -> str:
    candidates = [
        "genotype_file_path",
        "genotype_file",
        "file_path",
        "path",
    ]
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    raise ValueError(
        "Could not determine file column. Expected one of: "
        + ", ".join(candidates)
    )


def count_lines(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        return sum(1 for _ in handle)


def resolve_path(raw_value: str, assets_dir: Path, fallback_dir: Path) -> Path:
    candidate = Path(raw_value)
    if candidate.is_absolute():
        return candidate
    resolved = (assets_dir / candidate).resolve()
    if resolved.exists():
        return resolved
    return (fallback_dir / candidate).resolve()


def annotate_counts(input_path: Path, output_path: Path, assets_dir: Path) -> None:
    with input_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames: List[str] = reader.fieldnames or []
        if not fieldnames:
            raise ValueError("Input samplesheet has no header row")
        file_column = detect_file_column(fieldnames)
        base_dir = input_path.parent
        rows = []
        for row in reader:
            value = (row.get(file_column) or "").strip()
            if not value:
                row["line_count"] = "0"
                rows.append(row)
                continue
            candidate = resolve_path(value, assets_dir, base_dir)
            if candidate.exists():
                row["line_count"] = str(count_lines(candidate))
            else:
                row["line_count"] = "0"
            rows.append(row)

    output_fields = list(fieldnames)
    if "line_count" not in output_fields:
        output_fields.append("line_count")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).resolve()
    output_path = Path(args.output)
    assets_dir = Path(args.assets_dir).resolve()
    annotate_counts(input_path, output_path, assets_dir)


if __name__ == "__main__":
    main()
