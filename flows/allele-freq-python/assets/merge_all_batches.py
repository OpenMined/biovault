#!/usr/bin/env python3
import argparse
import os
from collections import defaultdict

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge dosage batch files into final outputs.")
    parser.add_argument("--batch-files", required=True, help="Space-delimited list of batch files")
    parser.add_argument("--matrix-tsv", required=True)
    parser.add_argument("--npz", required=True)
    parser.add_argument("--loci", required=True)
    parser.add_argument("--participants", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    batch_files = args.batch_files.split()

    participant_ids = set()
    loci_set = set()
    values = defaultdict(dict)
    rsid_map = {}

    for bf in batch_files:
        if not os.path.exists(bf):
            continue
        with open(bf, encoding="utf-8") as fh:
            next(fh, None)
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                locus = parts[0]
                rsid = parts[1]
                participant_id = parts[2]
                dosage = int(parts[3])
                participant_ids.add(participant_id)
                loci_set.add(locus)
                if rsid and rsid != ".":
                    if not rsid_map.get(locus):
                        rsid_map[locus] = rsid
                existing = values.get(locus, {}).get(participant_id)
                if existing is None or (existing == -1 and dosage != -1):
                    values.setdefault(locus, {})[participant_id] = dosage

    participants = sorted(participant_ids)
    sorted_loci = sorted(loci_set)
    n_loci = len(sorted_loci)
    n_participants = len(participants)

    idx = {locus: i for i, locus in enumerate(sorted_loci)}
    pidx = {pid: i for i, pid in enumerate(participants)}

    matrix = np.full((n_loci, n_participants), -1, dtype=np.int16)

    for locus, row in values.items():
        i = idx[locus]
        for pid, dosage in row.items():
            j = pidx[pid]
            matrix[i, j] = dosage

    with open(args.matrix_tsv, "w", encoding="utf-8") as out:
        out.write("locus_key\trsid\t" + "\t".join(participants) + "\n")
        for locus in sorted_loci:
            i = idx[locus]
            rsid = rsid_map.get(locus, "")
            out.write(f"{locus}\t{rsid}")
            for j in range(n_participants):
                out.write(f"\t{int(matrix[i, j])}")
            out.write("\n")

    np.savez_compressed(args.npz, dosage_mat=matrix, participants=participants, loci=sorted_loci)

    with open(args.participants, "w", encoding="utf-8") as out:
        out.write("participant_id\n")
        for pid in participants:
            out.write(f"{pid}\n")

    with open(args.loci, "w", encoding="utf-8") as out:
        out.write("#format=loci-v1\n")
        out.write("#genome=GRCh38\n")
        out.write(f"#n_loci={n_loci}\n")
        out.write("#key=CHROM-POS-REF-ALT\n")
        out.write("locus_key\n")
        for locus in sorted_loci:
            out.write(f"{locus}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
