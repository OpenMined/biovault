#!/usr/bin/env python3
"""
Validate TSV allele frequency outputs and compare client data.
"""
import argparse
import json
from pathlib import Path


def load_index(path: Path):
    data = json.loads(path.read_text())
    loci = [str(x) for x in data.get("loci", [])]
    rsids = [str(x) for x in data.get("rsids", [""] * len(loci))]
    if len(rsids) < len(loci):
        rsids.extend([""] * (len(loci) - len(rsids)))
    return loci, rsids


def load_tsv(path: Path):
    """Load TSV and return ordered loci, rsids, and ac values."""
    loci = []
    rsids = []
    ac_values = []
    with path.open() as f:
        header = f.readline().strip().split("\t")
        idx = {name: i for i, name in enumerate(header)}
        locus_col = idx.get("locus_key") if "locus_key" in idx else idx.get("locus")
        ac_col = idx.get("allele_count") if "allele_count" in idx else idx.get("ac")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            loci.append(parts[locus_col])
            rsids.append(parts[idx.get("rsid", 1)] if "rsid" in idx else "")
            ac_values.append(int(parts[ac_col]))
    return loci, rsids, ac_values


def check_client(label: str, base: Path):
    idx_file = base / "locus_index.json"
    tsv = base / "allele_freq.tsv"

    if not tsv.exists():
        print(f"[{label}] missing allele_freq.tsv")
        return None
    if not idx_file.exists():
        print(f"[{label}] missing locus_index.json")
        return None

    loci, rsids, ac = load_tsv(tsv)
    idx_loci, idx_rsids = load_index(idx_file)

    # Verify index matches TSV
    mismatches = 0
    for i, (t_locus, i_locus) in enumerate(zip(loci[:200], idx_loci[:200])):
        if t_locus != i_locus:
            mismatches += 1

    print(f"[{label}] loci={len(loci)} index_loci={len(idx_loci)} mismatches(sampled)={mismatches}")
    return {"ac": ac, "loci": loci, "rsids": rsids}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", required=True)
    parser.add_argument("--results-dir")
    args = parser.parse_args()

    workspace = Path(args.workspace)
    if args.results_dir:
        base = Path(args.results_dir)
    else:
        base = workspace / "results" / "flows" / "allele-freq-syqure"
    root = base / "gen_allele_freq"
    c1 = root / "client1_sandbox_local"
    c2 = root / "client2_sandbox_local"

    c1_data = check_client("client1", c1)
    c2_data = check_client("client2", c2)
    if not c1_data or not c2_data:
        print("missing client outputs or failed to decode npz")
        return 0

    apol1 = {"rs60910145", "rs71785313", "rs73885319"}
    thal = {"rs33985472", "rs33971634"}

    def report(data):
        rows = []
        for locus, rsid, ac in zip(data["loci"], data["rsids"], data["ac"]):
            if rsid in apol1 or rsid in thal:
                rows.append((rsid, locus, ac))
        return rows

    rows1 = report(c1_data)
    rows2 = report(c2_data)
    print("\nAPOL1/THAL rows (client1):")
    for rsid, locus, ac in rows1:
        print(f"  {rsid}\t{locus}\tac={ac}")
    print("\nAPOL1/THAL rows (client2):")
    for rsid, locus, ac in rows2:
        print(f"  {rsid}\t{locus}\tac={ac}")

    same = c1_data["ac"] == c2_data["ac"]
    print(f"\nAC arrays identical? {same}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
