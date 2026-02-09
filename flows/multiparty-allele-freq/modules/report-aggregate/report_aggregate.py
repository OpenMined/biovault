#!/usr/bin/env python3
"""
Generate aggregated allele frequency report from MPC aggregation results.

The input aggregated_counts is a concatenated vector: [ac..., an...]
We split it back into AC and AN halves to compute allele frequencies.
"""
import argparse
import json
import os
from pathlib import Path


def split_combined_counts(combined: list[int], loci_count: int) -> tuple[list[int], list[int]]:
    """
    Split a combined [ac..., an...] vector.
    If data is truncated (debug mode), fall back to AC-only with AN=zeros.
    """
    if loci_count <= 0:
        return [], []
    if len(combined) >= loci_count * 2:
        return combined[:loci_count], combined[loci_count : loci_count * 2]
    ac = combined[: min(len(combined), loci_count)]
    an = [0 for _ in range(len(ac))]
    return ac, an


def atomic_write_text(path: Path, content: str) -> None:
    """
    Write file content via temp file + atomic rename to avoid partial-file sync races.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with tmp_path.open("w", encoding="utf-8") as f:
        f.write(content)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp_path, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--union-index", required=True)
    parser.add_argument("--local-counts", required=True, help="Local aligned counts (for comparison)")
    parser.add_argument("--aggregated-counts", required=True, help="Combined AC+AN from MPC")
    parser.add_argument("--out-json")
    parser.add_argument("--out-tsv")
    parser.add_argument("--out-agg-tsv", help="Aggregated allele_freq.tsv format")
    parser.add_argument("--top", type=int, default=10)
    args = parser.parse_args()

    union_index = Path(args.union_index)
    local_counts_path = Path(args.local_counts)
    aggregated_counts_path = Path(args.aggregated_counts)

    for p in [union_index, local_counts_path, aggregated_counts_path]:
        if not p.exists():
            raise SystemExit(f"missing file: {p}")

    index = json.loads(union_index.read_text())
    loci = index.get("loci", [])
    rsids = index.get("rsids", [""] * len(loci))
    n = len(loci)

    # Load local counts (combined AC+AN)
    local_combined = json.loads(local_counts_path.read_text())
    local_ac, local_an = split_combined_counts(local_combined, n)

    # Load aggregated counts (combined AC+AN from MPC)
    agg_combined = json.loads(aggregated_counts_path.read_text())
    agg_ac, agg_an = split_combined_counts(agg_combined, n)

    # Ensure arrays match loci count
    if len(local_ac) != n or len(agg_ac) != n:
        print(f"WARNING: Array length mismatch. loci={n}, local_ac={len(local_ac)}, agg_ac={len(agg_ac)}")
        n = min(len(local_ac), len(agg_ac), n)
        local_ac = local_ac[:n]
        local_an = local_an[:n]
        agg_ac = agg_ac[:n]
        agg_an = agg_an[:n]
        if len(local_an) < n:
            local_an = local_an + [0 for _ in range(n - len(local_an))]
        if len(agg_an) < n:
            agg_an = agg_an + [0 for _ in range(n - len(agg_an))]
        loci = loci[:n]
        rsids = rsids[:n]

    total_agg_ac = sum(agg_ac)
    total_agg_an = sum(agg_an)
    total_local_ac = sum(local_ac)
    total_local_an = sum(local_an)
    nonzero_agg = sum(1 for x in agg_ac if x)
    max_agg = max(agg_ac) if agg_ac else 0

    print(f"aggregated_counts: {aggregated_counts_path}")
    print(f"loci={n} total_agg_ac={total_agg_ac} total_agg_an={total_agg_an}")
    print(f"total_local_ac={total_local_ac} total_local_an={total_local_an}")
    print(f"nonzero_agg={nonzero_agg} max_agg={max_agg}")

    top_idx = sorted(range(n), key=lambda i: agg_ac[i], reverse=True)[: min(args.top, n)]
    print("\nTop loci by aggregate count:")
    print("rank\tagg_ac\tagg_an\tagg_af\tlocal_ac\trsid\tlocus")
    for rank, i in enumerate(top_idx, start=1):
        ac = agg_ac[i]
        if ac == 0:
            continue
        an = agg_an[i]
        af = (ac / an) if an else 0.0
        local = local_ac[i]
        rs = rsids[i] if i < len(rsids) else ""
        locus = loci[i] if i < len(loci) else ""
        print(f"{rank}\t{ac}\t{an}\t{af:.6f}\t{local}\t{rs}\t{locus}")

    # Write detailed TSV report
    if args.out_tsv:
        out_path = Path(args.out_tsv)
        lines = ["locus\trsid\tagg_ac\tagg_an\tagg_af\tlocal_ac\tlocal_an\tlocal_af\n"]
        for i, (locus, rsid) in enumerate(zip(loci, rsids)):
            ac = agg_ac[i]
            an = agg_an[i]
            af = (ac / an) if an else 0.0
            lac = local_ac[i]
            lan = local_an[i]
            laf = (lac / lan) if lan else 0.0
            lines.append(f"{locus}\t{rsid}\t{ac}\t{an}\t{af:.6f}\t{lac}\t{lan}\t{laf:.6f}\n")
        atomic_write_text(out_path, "".join(lines))

    # Write aggregated allele_freq.tsv format (same format as input TSV but with aggregated values)
    if args.out_agg_tsv:
        out_path = Path(args.out_agg_tsv)
        lines = ["locus_key\trsid\tallele_count\tallele_number\tallele_freq\n"]
        for i, (locus, rsid) in enumerate(zip(loci, rsids)):
            ac = agg_ac[i]
            an = agg_an[i]
            af = (ac / an) if an else 0.0
            lines.append(f"{locus}\t{rsid}\t{ac}\t{an}\t{af:.6f}\n")
        atomic_write_text(out_path, "".join(lines))
        print(f"aggregated_allele_freq.tsv: {out_path}")

    # Write JSON summary
    if args.out_json:
        out_path = Path(args.out_json)
        top_rows = []
        for rank, i in enumerate(top_idx, start=1):
            ac = agg_ac[i]
            if ac == 0:
                continue
            an = agg_an[i]
            af = (ac / an) if an else 0.0
            local = local_ac[i]
            rs = rsids[i] if i < len(rsids) else ""
            locus = loci[i] if i < len(loci) else ""
            top_rows.append(
                {
                    "rank": rank,
                    "agg_ac": int(ac),
                    "agg_an": int(an),
                    "agg_af": af,
                    "local_ac": int(local),
                    "rsid": rs,
                    "locus": locus,
                }
            )
        payload = {
            "aggregated_counts_path": str(aggregated_counts_path),
            "loci": n,
            "total_agg_ac": int(total_agg_ac),
            "total_agg_an": int(total_agg_an),
            "total_local_ac": int(total_local_ac),
            "total_local_an": int(total_local_an),
            "nonzero_agg": int(nonzero_agg),
            "max_agg": int(max_agg),
            "top": top_rows,
        }
        atomic_write_text(out_path, json.dumps(payload, separators=(",", ":")) + "\n")

    print(f"report_json: {args.out_json}")
    print(f"report_tsv: {args.out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
