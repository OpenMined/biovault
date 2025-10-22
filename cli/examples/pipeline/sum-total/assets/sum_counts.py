#!/usr/bin/env python3
"""Sum the line_count column from a samplesheet and emit a text summary."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to the input samplesheet CSV")
    parser.add_argument("--output", required=True, help="Path to write the summary text file")
    return parser.parse_args()


def sum_counts(input_path: Path) -> int:
    with input_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "line_count" not in reader.fieldnames:
            raise ValueError("Input samplesheet must contain a 'line_count' column")
        total = 0
        for row in reader:
            value = row.get("line_count")
            try:
                total += int(value)
            except (TypeError, ValueError):
                continue
        return total


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).resolve()
    output_path = Path(args.output)
    total = sum_counts(input_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(f"total_line_count: {total}\n")


if __name__ == "__main__":
    main()
