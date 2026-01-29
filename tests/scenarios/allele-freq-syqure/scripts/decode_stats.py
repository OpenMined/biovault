#!/usr/bin/env python3
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
        / "05-secure_aggregate_ac"
        / "sum_ac.json"
    )
    sum_paths = [shared_glob] if shared_glob.exists() else []
    if not sum_paths:
        sum_paths = list(
            (workspace / "sandbox").glob(
                f"*/datasites/*/shared/flows/allele-freq-syqure/{run_id}/05-secure_aggregate_ac/sum_ac.json"
            )
        )
    if sum_paths:
        sum_path = sum_paths[0]
        agg = json.loads(sum_path.read_text())
    else:
        # Fall back to last aggregated_counts.json in results if MPC was skipped.
        fallback = list(base_results.glob("secure_aggregate_ac/*/aggregated_counts.json"))
        if not fallback:
            print("sum_ac.json not found for this run and no aggregated_counts.json fallback")
            return 0
        best = load_best_aggregate(fallback)
        if not best:
            print("aggregated_counts.json fallback found but failed to load")
            return 0
        sum_path, agg = best

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

    c1 = load_client(
        base_results
        / "align_counts"
        / "client1_sandbox_local"
        / "aligned_ac.json"
    )
    c2 = load_client(
        base_results
        / "align_counts"
        / "client2_sandbox_local"
        / "aligned_ac.json"
    )

    n = len(agg)
    print(f"sum_ac.json: {sum_path}")
    print(
        f"loci: {len(loci)} agg_len: {n} c1_len: {len(c1) if c1 else 'n/a'} c2_len: {len(c2) if c2 else 'n/a'}"
    )

    total_agg = sum(agg)
    nz = sum(1 for x in agg if x)
    max_val = max(agg) if agg else 0
    print(f"total_agg={total_agg} nonzero={nz} max={max_val}")
    if c1:
        total_c1 = sum(c1)
        print(
            f"total_c1={total_c1} share_vs_agg={total_c1 / total_agg * 100:.2f}%"
            if total_agg
            else "total_c1=0"
        )
    if c2:
        total_c2 = sum(c2)
        print(
            f"total_c2={total_c2} share_vs_agg={total_c2 / total_agg * 100:.2f}%"
            if total_agg
            else "total_c2=0"
        )

    top_n = 10
    top_idx = sorted(range(n), key=lambda i: agg[i], reverse=True)[:top_n]
    print("\nTop loci by aggregate count:")
    print("rank\tagg\tc1\tc2\t%c1\t%c2\trsid\tlocus")
    for rank, i in enumerate(top_idx, start=1):
        a = agg[i]
        if a == 0:
            continue
        v1 = c1[i] if c1 else None
        v2 = c2[i] if c2 else None
        p1 = (v1 / a * 100.0) if (v1 is not None and a) else float("nan")
        p2 = (v2 / a * 100.0) if (v2 is not None and a) else float("nan")
        rs = rsids[i] if i < len(rsids) else ""
        locus = loci[i] if i < len(loci) else ""
        print(
            f"{rank}\t{a}\t{v1 if v1 is not None else '-'}\t{v2 if v2 is not None else '-'}\t{p1:.2f}\t{p2:.2f}\t{rs}\t{locus}"
        )

    if c1 and c2:
        diffs = []
        for i, a in enumerate(agg):
            if a:
                diffs.append(abs(c1[i] - c2[i]) / a * 100.0)
        diffs.sort(reverse=True)
        if diffs:
            print(f"\nMax % diff between clients: {diffs[0]:.2f}%")
            print(f"Mean % diff between clients: {sum(diffs)/len(diffs):.2f}%")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
