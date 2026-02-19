#!/usr/bin/env python3
"""
Generate aggregated allele frequency report from MPC aggregation results.

Supports two aggregated_counts formats:
  - New (HE+Chebyshev): [af_0, ..., af_n] — pre-computed allele frequencies (floats, length n)
  - Legacy (SMPC):      [ac_0,...,ac_n, an_0,...,an_n] — combined counts (ints, length 2n)

Detection: if len(data) == n, treat as AF floats; if len(data) == 2*n, split into AC/AN.
"""
import argparse
import json
import os
from pathlib import Path


def split_combined_counts(combined: list, loci_count: int) -> tuple[list[int], list[int]]:
    if loci_count <= 0:
        return [], []
    if len(combined) >= loci_count * 2:
        return combined[:loci_count], combined[loci_count : loci_count * 2]
    ac = combined[: min(len(combined), loci_count)]
    an = [0 for _ in range(len(ac))]
    return ac, an


def atomic_write_text(path: Path, content: str) -> None:
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
    parser.add_argument("--aggregated-counts", required=True, help="Aggregated output from MPC")
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

    # Load local counts (always combined AC+AN from align_counts)
    local_combined = json.loads(local_counts_path.read_text())
    local_ac, local_an = split_combined_counts(local_combined, n)

    # Load aggregated data — detect format
    agg_raw = json.loads(aggregated_counts_path.read_text())

    if len(agg_raw) == n:
        # New format: pre-computed allele frequencies from HE+Chebyshev
        agg_af_raw = [float(x) for x in agg_raw]
        out_of_range = sum(1 for x in agg_af_raw if x < 0.0 or x > 1.0)
        agg_af = [max(0.0, min(1.0, x)) for x in agg_af_raw]
        if out_of_range > 0:
            print(f"Clamped {out_of_range}/{n} AF values to [0, 1]")
        agg_ac = None
        agg_an = None
        mode = "af_direct"
        print(f"aggregated format: AF direct (HE+Chebyshev), n={n}")
    elif len(agg_raw) >= n * 2:
        # Legacy format: combined [ac..., an...] ints
        agg_ac_raw, agg_an_raw = split_combined_counts(agg_raw, n)
        agg_ac = [int(x) for x in agg_ac_raw]
        agg_an = [int(x) for x in agg_an_raw]
        agg_af = [(ac / an) if an else 0.0 for ac, an in zip(agg_ac, agg_an)]
        mode = "ac_an_split"
        print(f"aggregated format: AC+AN split (legacy SMPC), n={n}")
    else:
        # Fallback: treat as AF if shorter than 2*n
        agg_af = [float(x) for x in agg_raw] + [0.0] * max(0, n - len(agg_raw))
        agg_ac = None
        agg_an = None
        mode = "af_direct"
        print(f"WARNING: aggregated array length {len(agg_raw)} != n={n} or 2n={2*n}; treating as AF")

    # Trim to common length
    if len(local_ac) != n or len(agg_af) != n:
        print(f"WARNING: Array length mismatch. loci={n}, local_ac={len(local_ac)}, agg_af={len(agg_af)}")
        n = min(len(local_ac), len(agg_af), n)
        local_ac = local_ac[:n]
        local_an = local_an[:n]
        agg_af = agg_af[:n]
        if agg_ac is not None:
            agg_ac = agg_ac[:n]
            agg_an = agg_an[:n]
        loci = loci[:n]
        rsids = rsids[:n]

    total_local_ac = sum(local_ac)
    total_local_an = sum(local_an)
    nonzero_af = sum(1 for x in agg_af if x > 0.0)
    max_af = max(agg_af) if agg_af else 0.0

    print(f"aggregated_counts: {aggregated_counts_path}")
    print(f"loci={n} mode={mode}")
    if agg_ac is not None:
        print(f"total_agg_ac={sum(agg_ac)} total_agg_an={sum(agg_an)}")
    print(f"total_local_ac={total_local_ac} total_local_an={total_local_an}")
    print(f"nonzero_af={nonzero_af} max_af={max_af:.6f}")

    top_idx = sorted(range(n), key=lambda i: agg_af[i], reverse=True)[: min(args.top, n)]
    print("\nTop loci by aggregated allele frequency:")
    print("rank\tagg_af\tlocal_ac\tlocal_an\tlocal_af\trsid\tlocus")
    for rank, i in enumerate(top_idx, start=1):
        af = agg_af[i]
        if af == 0.0:
            continue
        lac = local_ac[i]
        lan = local_an[i]
        laf = (lac / lan) if lan else 0.0
        rs = rsids[i] if i < len(rsids) else ""
        locus = loci[i] if i < len(loci) else ""
        print(f"{rank}\t{af:.6f}\t{lac}\t{lan}\t{laf:.6f}\t{rs}\t{locus}")

    # Write detailed TSV report
    if args.out_tsv:
        out_path = Path(args.out_tsv)
        lines = ["locus\trsid\tagg_af\tlocal_ac\tlocal_an\tlocal_af\n"]
        for i, (locus, rsid) in enumerate(zip(loci, rsids)):
            af = agg_af[i]
            lac = local_ac[i]
            lan = local_an[i]
            laf = (lac / lan) if lan else 0.0
            lines.append(f"{locus}\t{rsid}\t{af:.6f}\t{lac}\t{lan}\t{laf:.6f}\n")
        atomic_write_text(out_path, "".join(lines))

    # Write aggregated allele_freq.tsv format
    if args.out_agg_tsv:
        out_path = Path(args.out_agg_tsv)
        lines = ["locus_key\trsid\tallele_freq\n"]
        for i, (locus, rsid) in enumerate(zip(loci, rsids)):
            af = agg_af[i]
            lines.append(f"{locus}\t{rsid}\t{af:.6f}\n")
        atomic_write_text(out_path, "".join(lines))
        print(f"aggregated_allele_freq.tsv: {out_path}")

    # Write JSON summary
    if args.out_json:
        out_path = Path(args.out_json)
        top_rows = []
        for rank, i in enumerate(top_idx, start=1):
            af = agg_af[i]
            if af == 0.0:
                continue
            lac = local_ac[i]
            lan = local_an[i]
            laf = (lac / lan) if lan else 0.0
            rs = rsids[i] if i < len(rsids) else ""
            locus = loci[i] if i < len(loci) else ""
            row = {
                "rank": rank,
                "agg_af": round(af, 6),
                "local_ac": int(lac),
                "local_an": int(lan),
                "local_af": round(laf, 6),
                "rsid": rs,
                "locus": locus,
            }
            if agg_ac is not None:
                row["agg_ac"] = int(agg_ac[i])
                row["agg_an"] = int(agg_an[i])
            top_rows.append(row)
        payload = {
            "aggregated_counts_path": str(aggregated_counts_path),
            "mode": mode,
            "loci": n,
            "total_local_ac": int(total_local_ac),
            "total_local_an": int(total_local_an),
            "nonzero_af": int(nonzero_af),
            "max_af": round(max_af, 6),
            "top": top_rows,
        }
        if agg_ac is not None:
            payload["total_agg_ac"] = int(sum(agg_ac))
            payload["total_agg_an"] = int(sum(agg_an))
        atomic_write_text(out_path, json.dumps(payload, separators=(",", ":")) + "\n")

    print(f"report_json: {args.out_json}")
    print(f"report_tsv: {args.out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
