#!/usr/bin/env python3
"""
Trusted aggregation of allele frequency TSVs from multiple clients.

Reads a manifest of client TSV paths, loads each, computes the union locus
index, aligns all datasets, sums AC/AN, and outputs aggregated results.
"""
import argparse
import json
import os
from pathlib import Path

from allele_freq_utils import AlleleFreqData, LocusIndex, aggregate


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
    parser.add_argument("--manifest", required=True, help="Manifest file listing client TSV paths")
    parser.add_argument("--out-tsv", required=True, help="Output aggregated allele_freq.tsv")
    parser.add_argument("--out-json", required=True, help="Output report.json")
    parser.add_argument("--out-index", required=True, help="Output union_locus_index.json")
    parser.add_argument("--top", type=int, default=10)
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    tsv_paths = []
    for line in manifest_path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        # Manifest format: "email\tpath" — extract path (second column)
        parts = line.split("\t")
        tsv_paths.append(parts[-1].strip())
    print(f"manifest: {manifest_path} ({len(tsv_paths)} client TSVs)")

    datasets = []
    for tsv_path in tsv_paths:
        p = Path(tsv_path)
        if not p.exists():
            print(f"WARNING: missing client TSV: {p}")
            continue
        data = AlleleFreqData.from_tsv(str(p))
        print(f"  loaded {p.name}: {len(data)} loci, total_ac={int(data.ac.sum())}, total_an={int(data.an.sum())}")
        datasets.append(data)

    if len(datasets) < 2:
        raise SystemExit(f"Need at least 2 client datasets, got {len(datasets)}")

    agg = aggregate(*datasets)
    n = len(agg)
    af = agg.af
    print(f"aggregated: {n} loci, total_ac={int(agg.ac.sum())}, total_an={int(agg.an.sum())}")

    # Save union locus index
    out_index = Path(args.out_index)
    agg.index.save(str(out_index))
    print(f"union_index: {out_index} ({n} loci)")

    # Save aggregated allele_freq TSV
    out_tsv = Path(args.out_tsv)
    agg.save_tsv(str(out_tsv))
    print(f"aggregated_allele_freq.tsv: {out_tsv} ({n} loci)")

    # Build report JSON
    nonzero_af = int(sum(1 for x in af if x > 0.0))
    max_af = float(max(af)) if len(af) > 0 else 0.0
    top_idx = sorted(range(n), key=lambda i: af[i], reverse=True)[:min(args.top, n)]

    print(f"\nnonzero_af={nonzero_af} max_af={max_af:.6f}")
    print("\nTop loci by aggregated allele frequency:")
    print("rank\taf\tac\tan\tn_samples\trsid\tlocus")

    top_rows = []
    for rank, i in enumerate(top_idx, start=1):
        if af[i] == 0.0:
            continue
        row = {
            "rank": rank,
            "af": round(float(af[i]), 6),
            "ac": int(agg.ac[i]),
            "an": int(agg.an[i]),
            "n_samples": int(agg.n_samples[i]),
            "rsid": str(agg.index.rsids[i]),
            "locus": str(agg.index.loci[i]),
        }
        top_rows.append(row)
        print(f"{rank}\t{af[i]:.6f}\t{agg.ac[i]}\t{agg.an[i]}\t{agg.n_samples[i]}\t{agg.index.rsids[i]}\t{agg.index.loci[i]}")

    payload = {
        "mode": "trusted_plaintext",
        "n_clients": len(datasets),
        "loci": n,
        "total_ac": int(agg.ac.sum()),
        "total_an": int(agg.an.sum()),
        "nonzero_af": nonzero_af,
        "max_af": round(max_af, 6),
        "top": top_rows,
    }

    out_json = Path(args.out_json)
    atomic_write_text(out_json, json.dumps(payload, indent=2) + "\n")
    print(f"report_json: {out_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
