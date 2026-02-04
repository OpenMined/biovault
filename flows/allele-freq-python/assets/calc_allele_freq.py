#!/usr/bin/env python3
import argparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Calculate allele frequencies from dosage matrix.")
    parser.add_argument("--matrix", required=True, help="Dosage matrix TSV")
    parser.add_argument("--output", required=True, help="Output allele frequency TSV")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with open(args.matrix, encoding="utf-8") as fh:
        header = next(fh).strip().split("\t")
        if len(header) < 3 or header[0] != "locus_key" or header[1] != "rsid":
            raise SystemExit("Unexpected matrix header; expected: locus_key\\trsid\\t<participants...>")
        pids = header[2:]
        rows = [line.strip().split("\t") for line in fh if line.strip()]

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("locus_key\tallele_count\tallele_number\tnum_homo\tallele_freq\trsid\n")

        for row in rows:
            if len(row) < 2:
                continue
            locus = row[0]
            rsid = row[1]
            vals = [int(v) for v in row[2:]]
            observed = [v for v in vals if v != -1]
            n_obs = len(observed)
            allele_number = 2 * n_obs
            allele_count = sum(observed)
            num_homo = sum(1 for v in observed if v == 2)
            af = allele_count / allele_number if allele_number > 0 else 0.0
            out.write(f"{locus}\t{allele_count}\t{allele_number}\t{num_homo}\t{af:.6f}\t{rsid}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
