#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--union-index", required=True)
    parser.add_argument("--aligned-ac", required=True)
    parser.add_argument("--aggregated", required=True)
    parser.add_argument("--out-json")
    parser.add_argument("--out-tsv")
    parser.add_argument("--top", type=int, default=10)
    args = parser.parse_args()

    union_index = Path(args.union_index)
    aligned_ac = Path(args.aligned_ac)
    aggregated = Path(args.aggregated)
    if not union_index.exists() or not aligned_ac.exists() or not aggregated.exists():
        raise SystemExit("missing union_index, aligned_ac, or aggregated_counts")

    index = json.loads(union_index.read_text())
    loci = index.get("loci", [])
    rsids = index.get("rsids", [""] * len(loci))
    local_ac = json.loads(aligned_ac.read_text())
    agg_ac = json.loads(aggregated.read_text())

    n = min(len(local_ac), len(agg_ac), len(loci))
    local_ac = local_ac[:n]
    agg_ac = agg_ac[:n]
    loci = loci[:n]
    rsids = rsids[:n]

    total_agg = sum(agg_ac)
    total_local = sum(local_ac)
    nonzero_agg = sum(1 for x in agg_ac if x)
    max_agg = max(agg_ac) if agg_ac else 0

    print(f"aggregate: {aggregated}")
    print(f"loci={n} total_agg={total_agg} total_local={total_local} nonzero_agg={nonzero_agg} max_agg={max_agg}")

    top_idx = sorted(range(n), key=lambda i: agg_ac[i], reverse=True)[: min(args.top, n)]
    print("\nTop loci by aggregate count:")
    print("rank\tagg\tlocal\t%local\trsid\tlocus")
    for rank, i in enumerate(top_idx, start=1):
        agg = agg_ac[i]
        if agg == 0:
            continue
        local = local_ac[i]
        pct = (local / agg * 100.0) if agg else 0.0
        rs = rsids[i] if i < len(rsids) else ""
        locus = loci[i] if i < len(loci) else ""
        print(f"{rank}\t{agg}\t{local}\t{pct:.2f}\t{rs}\t{locus}")

    if args.out_tsv:
        out_path = Path(args.out_tsv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8") as f:
            f.write("locus\trsid\tlocal_ac\tagg_ac\tpct_local_of_agg\tdiff\n")
            for locus, rsid, local, agg in zip(loci, rsids, local_ac, agg_ac):
                if agg:
                    pct = local / agg * 100.0
                else:
                    pct = 0.0
                diff = agg - local
                f.write(f"{locus}\t{rsid}\t{local}\t{agg}\t{pct:.6f}\t{diff}\n")

    if args.out_json:
        out_path = Path(args.out_json)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        top_rows = []
        for rank, i in enumerate(top_idx, start=1):
            agg = agg_ac[i]
            if agg == 0:
                continue
            local = local_ac[i]
            pct = (local / agg * 100.0) if agg else 0.0
            rs = rsids[i] if i < len(rsids) else ""
            locus = loci[i] if i < len(loci) else ""
            top_rows.append(
                {
                    "rank": rank,
                    "agg": int(agg),
                    "local": int(local),
                    "pct_local": pct,
                    "rsid": rs,
                    "locus": locus,
                }
            )
        payload = {
            "aggregate_path": str(aggregated),
            "aligned_ac_path": str(aligned_ac),
            "loci": n,
            "total_agg": int(total_agg),
            "total_local": int(total_local),
            "nonzero_agg": int(nonzero_agg),
            "max_agg": int(max_agg),
            "top": top_rows,
        }
        out_path.write_text(json.dumps(payload, separators=(",", ":")) + "\n", encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
