#!/usr/bin/env python3
"""
Verify that the MPC aggregation produces correct results.

This script:
1. Loads each client's allele_freq.tsv
2. Aligns them to the union index
3. Sums the AC and AN columns
4. Compares to the aggregated_allele_freq.tsv values
5. Asserts they match exactly

Usage:
    python verify_aggregation.py --results-dir <path> --workspace <path>
"""
import argparse
import json
import sys
from pathlib import Path


def load_tsv(path: Path) -> dict:
    """Load a TSV file into a dict keyed by locus_key."""
    data = {}
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        # Handle both old and new column names
        locus_col = col_idx.get("locus_key") if "locus_key" in col_idx else col_idx.get("locus")
        ac_col = col_idx.get("allele_count") if "allele_count" in col_idx else col_idx.get("ac", col_idx.get("agg_ac"))
        an_col = col_idx.get("allele_number") if "allele_number" in col_idx else col_idx.get("an", col_idx.get("agg_an"))
        rsid_col = col_idx.get("rsid")

        if locus_col is None or ac_col is None or an_col is None:
            raise ValueError(f"Missing required columns in {path}. Header: {header}")

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= max(locus_col, ac_col, an_col):
                continue
            locus = parts[locus_col]
            ac = int(parts[ac_col])
            an = int(parts[an_col])
            rsid = parts[rsid_col] if rsid_col is not None and rsid_col < len(parts) else ""
            data[locus] = {"ac": ac, "an": an, "rsid": rsid}
    return data


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify MPC aggregation correctness")
    parser.add_argument("--results-dir", required=True, help="Flow results directory")
    parser.add_argument("--workspace", required=True, help="Workspace root")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    workspace = Path(args.workspace)

    # Find client allele_freq.tsv files
    client1_tsv = results_dir / "gen_allele_freq" / "client1_sandbox_local" / "allele_freq.tsv"
    client2_tsv = results_dir / "gen_allele_freq" / "client2_sandbox_local" / "allele_freq.tsv"

    # Find aggregated allele_freq.tsv (from either client's report_aggregate output)
    agg_tsv = results_dir / "report_aggregate" / "client1_sandbox_local" / "aggregated_allele_freq.tsv"

    # Load union index for locus list
    union_index_path = results_dir / "build_master" / "aggregator_sandbox_local" / "union_locus_index.json"

    print(f"Client1 TSV: {client1_tsv}")
    print(f"Client2 TSV: {client2_tsv}")
    print(f"Aggregated TSV: {agg_tsv}")
    print(f"Union index: {union_index_path}")

    for p in [client1_tsv, client2_tsv, agg_tsv, union_index_path]:
        if not p.exists():
            print(f"ERROR: Missing file: {p}")
            return 1

    # Load data
    print("\nLoading client data...")
    client1_data = load_tsv(client1_tsv)
    client2_data = load_tsv(client2_tsv)
    agg_data = load_tsv(agg_tsv)

    with union_index_path.open("r") as f:
        union_index = json.load(f)

    loci = union_index.get("loci", [])
    print(f"Union index has {len(loci)} loci")
    print(f"Client1 has {len(client1_data)} loci")
    print(f"Client2 has {len(client2_data)} loci")
    print(f"Aggregated has {len(agg_data)} loci")

    # Verify aggregation for each locus
    mismatches = []
    total_expected_ac = 0
    total_expected_an = 0
    total_actual_ac = 0
    total_actual_an = 0

    for locus in loci:
        c1 = client1_data.get(locus, {"ac": 0, "an": 0})
        c2 = client2_data.get(locus, {"ac": 0, "an": 0})
        agg = agg_data.get(locus, {"ac": -1, "an": -1})

        expected_ac = c1["ac"] + c2["ac"]
        expected_an = c1["an"] + c2["an"]
        actual_ac = agg["ac"]
        actual_an = agg["an"]

        total_expected_ac += expected_ac
        total_expected_an += expected_an
        total_actual_ac += actual_ac
        total_actual_an += actual_an

        if expected_ac != actual_ac or expected_an != actual_an:
            mismatches.append({
                "locus": locus,
                "c1_ac": c1["ac"], "c1_an": c1["an"],
                "c2_ac": c2["ac"], "c2_an": c2["an"],
                "expected_ac": expected_ac, "expected_an": expected_an,
                "actual_ac": actual_ac, "actual_an": actual_an,
            })

    print(f"\nTotal expected AC: {total_expected_ac}")
    print(f"Total actual AC:   {total_actual_ac}")
    print(f"Total expected AN: {total_expected_an}")
    print(f"Total actual AN:   {total_actual_an}")

    if mismatches:
        print(f"\n❌ VERIFICATION FAILED: {len(mismatches)} mismatches found!")
        print("\nFirst 10 mismatches:")
        for m in mismatches[:10]:
            print(f"  {m['locus']}: c1=({m['c1_ac']},{m['c1_an']}) c2=({m['c2_ac']},{m['c2_an']}) "
                  f"expected=({m['expected_ac']},{m['expected_an']}) actual=({m['actual_ac']},{m['actual_an']})")
        return 1

    print(f"\n✅ VERIFICATION PASSED: All {len(loci)} loci match!")
    print("   client1.allele_freq.tsv + client2.allele_freq.tsv == aggregated_allele_freq.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
