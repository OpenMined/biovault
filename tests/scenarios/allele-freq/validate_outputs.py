#!/usr/bin/env python3
import sys
from pathlib import Path


REQUIRED_FILES = [
    "dosage_matrix.tsv",
    "dosage_matrix.npz",
    "locus_index.txt",
    "participants.txt",
    "allele_freq.tsv",
    "vcf_conversion_results.tsv",
]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: validate_outputs.py <results_dir>", file=sys.stderr)
        return 2

    results_dir = Path(sys.argv[1]).expanduser().resolve()
    if not results_dir.exists():
        print(f"results dir not found: {results_dir}", file=sys.stderr)
        return 2

    missing = []
    empty = []
    for name in REQUIRED_FILES:
        path = results_dir / name
        if not path.exists():
            missing.append(name)
            continue
        if path.is_file() and path.stat().st_size == 0:
            empty.append(name)

    if missing:
        candidate_dirs = {}
        for name in missing:
            for found in results_dir.rglob(name):
                candidate_dirs.setdefault(found.parent, set()).add(name)
        for directory, found_names in candidate_dirs.items():
            if len(found_names) == len(missing):
                missing = []
                for name in REQUIRED_FILES:
                    path = directory / name
                    if not path.exists():
                        missing.append(name)
                        continue
                    if path.is_file() and path.stat().st_size == 0:
                        empty.append(name)
                if not missing:
                    print(f"ok: outputs found under {directory}")
                    return 0

    if missing:
        print(f"missing outputs: {', '.join(missing)}", file=sys.stderr)
        return 1
    if empty:
        print(f"empty outputs: {', '.join(empty)}", file=sys.stderr)
        return 1

    print("ok: all required outputs present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
