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
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Base directory where referenced files must reside",
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


def resolve_candidate(value: str, data_dir: Path) -> Path | None:
    raw = Path(value)
    if raw.is_absolute():
        return raw if raw.exists() else None
    candidate = (data_dir / raw).resolve()
    return candidate if candidate.exists() else None


def write_filtered(
    input_path: Path,
    output_path: Path,
    extension: str,
    data_dir: Path,
) -> None:
    with input_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames: List[str] = reader.fieldnames or []
        if not fieldnames:
            raise ValueError("Input samplesheet has no header row")
        file_column = detect_file_column(fieldnames)

        rows: List[dict] = []
        for row in reader:
            value = (row.get(file_column) or "").strip()
            if not value:
                continue
            candidate = resolve_candidate(value, data_dir)
            if candidate is None:
                continue
            if candidate.suffix.lower() != extension.lower():
                continue
            # Preserve original row but normalise to relative path when possible
            updated = dict(row)
            if not Path(value).is_absolute():
                # Keep the original relative entry
                updated[file_column] = value
            else:
                try:
                    updated[file_column] = candidate.relative_to(data_dir).as_posix()
                except ValueError:
                    updated[file_column] = candidate.as_posix()
            rows.append(updated)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).resolve()
    output_path = Path(args.output)
    data_dir = Path(args.data_dir).resolve()
    write_filtered(input_path, output_path, args.extension, data_dir)


if __name__ == "__main__":
    main()
