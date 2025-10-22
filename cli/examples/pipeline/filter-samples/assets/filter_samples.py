#!/usr/bin/env python3
"""Filter samplesheet rows by file extension."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to the input samplesheet CSV")
    parser.add_argument("--output", required=True, help="Path to write the filtered CSV")
    parser.add_argument(
        "--extension",
        default=".txt",
        help="Keep rows whose file path ends with this extension (case-insensitive)",
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


def keep_row(row: dict, column: str, extension: str, base_dir: Path) -> tuple[bool, Path | None]:
    value = (row.get(column) or "").strip()
    if not value:
        return False, None
    candidate = Path(value)
    if not candidate.is_absolute():
        candidate = (base_dir / candidate).resolve()
    suffix = candidate.suffix.lower()
    if suffix == extension.lower() and candidate.exists():
        return True, candidate
    return False, None


def write_filtered(input_path: Path, output_path: Path, extension: str) -> None:
    with input_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames: List[str] = reader.fieldnames or []
        if not fieldnames:
            raise ValueError("Input samplesheet has no header row")
        file_column = detect_file_column(fieldnames)
        base_dir = input_path.parent

        rows: List[dict] = []
        for row in reader:
            keep, resolved = keep_row(row, file_column, extension, base_dir)
            if keep and resolved:
                updated_row = dict(row)
                updated_row[file_column] = str(resolved)
                rows.append(updated_row)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).resolve()
    output_path = Path(args.output)
    write_filtered(input_path, output_path, args.extension)


if __name__ == "__main__":
    main()
