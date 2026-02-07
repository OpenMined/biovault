#!/usr/bin/env python3
"""
Decode and display statistics from the MPC aggregation results.

Updated to handle combined AC+AN vectors (single MPC aggregation).
"""
import argparse
import json
from pathlib import Path


def load_client(path: Path):
    if not path.exists():
        return None
    return json.loads(path.read_text())


def load_best_aggregate(paths):
    best = None
    best_total = -1
    for path in paths:
        try:
            data = json.loads(path.read_text())
        except Exception:
            continue
        total = sum(data)
        if total > best_total:
            best_total = total
            best = (path, data)
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--results-dir")
    args = parser.parse_args()

    workspace = Path(args.workspace)
    run_id = args.run_id
    if args.results_dir:
        base_results = Path(args.results_dir)
    else:
        base_results = workspace / "results" / "flows" / "allele-freq-syqure"

    # Look for aggregated_counts.json in results (combined AC+AN vector)
    agg_paths = list(base_results.glob("secure_aggregate/*/aggregated_counts.json"))
    if not agg_paths:
        # Try shared path
        shared_glob = (
            workspace
            / "sandbox"
            / "client2@sandbox.local"
            / "datasites"
            / "client1@sandbox.local"
            / "shared"
            / "flows"
            / "allele-freq-syqure"
            / run_id
            / "*-secure_aggregate"
            / "aggregated_counts.json"
        )
        agg_paths = list(shared_glob.parent.parent.glob("*-secure_aggregate/aggregated_counts.json"))

    if not agg_paths:
        print("aggregated_counts.json not found")
        return 0

    best = load_best_aggregate(agg_paths)
    if not best:
        print("Failed to load aggregated_counts.json")
        return 0
    sum_path, agg_combined = best

    index_path = (
        base_results
        / "build_master"
        / "aggregator_sandbox_local"
        / "union_locus_index.json"
    )
    if not index_path.exists():
        raise SystemExit(f"Missing union_locus_index.json at {index_path}")
    index = json.loads(index_path.read_text())
    loci = index.get("loci", [])
    rsids = index.get("rsids", [""] * len(loci))
    n = len(loci)

    # Split combined vector into AC and AN halves
    agg_ac = agg_combined[:n]
    agg_an = agg_combined[n:] if len(agg_combined) > n else [0] * n

    # Load client aligned_counts.json (combined AC+AN vectors)
    c1_combined = load_client(
        base_results
        / "align_counts"
        / "client1_sandbox_local"
        / "aligned_counts.json"
    )
    c2_combined = load_client(
        base_results
        / "align_counts"
        / "client2_sandbox_local"
        / "aligned_counts.json"
    )

    # Split client vectors
    c1_ac = c1_combined[:n] if c1_combined else None
    c2_ac = c2_combined[:n] if c2_combined else None

    print(f"aggregated_counts.json: {sum_path}")
    print(
        f"loci: {len(loci)} agg_len: {len(agg_combined)} c1_len: {len(c1_combined) if c1_combined else 'n/a'} c2_len: {len(c2_combined) if c2_combined else 'n/a'}"
    )

    total_agg_ac = sum(agg_ac)
    total_agg_an = sum(agg_an)
    nz = sum(1 for x in agg_ac if x)
    max_val = max(agg_ac) if agg_ac else 0
    print(f"total_agg_ac={total_agg_ac} total_agg_an={total_agg_an} nonzero={nz} max={max_val}")

    if c1_ac:
        total_c1 = sum(c1_ac)
        print(
            f"total_c1={total_c1} share_vs_agg={total_c1 / total_agg_ac * 100:.2f}%"
            if total_agg_ac
            else "total_c1=0"
        )
    if c2_ac:
        total_c2 = sum(c2_ac)
        print(
            f"total_c2={total_c2} share_vs_agg={total_c2 / total_agg_ac * 100:.2f}%"
            if total_agg_ac
            else "total_c2=0"
        )

    top_n = 10
    top_idx = sorted(range(n), key=lambda i: agg_ac[i], reverse=True)[:top_n]
    print("\nTop loci by aggregate count:")
    print("rank\tagg_ac\tagg_an\tagg_af\tc1\tc2\t%c1\t%c2\trsid\tlocus")
    for rank, i in enumerate(top_idx, start=1):
        ac = agg_ac[i]
        if ac == 0:
            continue
        an = agg_an[i]
        af = (ac / an) if an else 0.0
        v1 = c1_ac[i] if c1_ac else None
        v2 = c2_ac[i] if c2_ac else None
        p1 = (v1 / ac * 100.0) if (v1 is not None and ac) else float("nan")
        p2 = (v2 / ac * 100.0) if (v2 is not None and ac) else float("nan")
        rs = rsids[i] if i < len(rsids) else ""
        locus = loci[i] if i < len(loci) else ""
        print(
            f"{rank}\t{ac}\t{an}\t{af:.4f}\t{v1 if v1 is not None else '-'}\t{v2 if v2 is not None else '-'}\t{p1:.2f}\t{p2:.2f}\t{rs}\t{locus}"
        )

    if c1_ac and c2_ac:
        diffs = []
        for i, ac in enumerate(agg_ac):
            if ac:
                diffs.append(abs(c1_ac[i] - c2_ac[i]) / ac * 100.0)
        diffs.sort(reverse=True)
        if diffs:
            print(f"\nMax % diff between clients: {diffs[0]:.2f}%")
            print(f"Mean % diff between clients: {sum(diffs)/len(diffs):.2f}%")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
