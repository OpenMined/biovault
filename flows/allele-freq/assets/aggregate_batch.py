#!/usr/bin/env python3
import argparse
import os


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Concatenate per-participant dosage files into a batch.")
    parser.add_argument("--batch-id", required=True)
    parser.add_argument("--counts", required=True, help="Space-delimited list of count files")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    count_files = args.counts.split()

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("locus_key\trsid\tparticipant_id\tdosage\n")
        for cf in count_files:
            if not os.path.exists(cf):
                continue
            with open(cf, encoding="utf-8") as fh:
                next(fh, None)
                for line in fh:
                    out.write(line)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
