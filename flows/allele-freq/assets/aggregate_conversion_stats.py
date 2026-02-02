#!/usr/bin/env python3
import argparse
import os


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate per-participant VCF conversion stats.")
    parser.add_argument("--stats-files", required=True, help="Space-delimited list of stats files")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    stats_files = args.stats_files.split()

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("participant_id\tfilename\tvcf_rows\tmissing_count\tinferred_count\tstatus\n")
        for sf in stats_files:
            if not os.path.exists(sf):
                continue
            with open(sf, encoding="utf-8") as fh:
                next(fh, None)
                for line in fh:
                    out.write(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
