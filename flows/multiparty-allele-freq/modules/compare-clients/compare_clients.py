#!/usr/bin/env python3
import json
import os
from pathlib import Path


APOL1_RSIDS = {"rs60910145", "rs71785313", "rs73885319"}
THAL_RSIDS = {"rs33985472", "rs33971634"}


def load_manifest(path: Path):
    entries = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        ds, p = parts[0].strip(), parts[1].strip()
        if ds and p:
            entries.append((ds, Path(p)))
    return entries


def load_index(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    loci = data.get("loci", [])
    rsids = data.get("rsids", [""] * len(loci))
    if len(rsids) < len(loci):
        rsids.extend([""] * (len(loci) - len(rsids)))
    return loci, rsids


def load_ac(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def top_rows(ac, loci, rsids, n=10):
    pairs = sorted(range(len(ac)), key=lambda i: ac[i], reverse=True)[:n]
    rows = []
    for i in pairs:
        rows.append((ac[i], rsids[i] if i < len(rsids) else "", loci[i] if i < len(loci) else ""))
    return rows


def overlay_rows(ac, loci, rsids):
    rows = []
    for i, rsid in enumerate(rsids):
        if rsid in APOL1_RSIDS or rsid in THAL_RSIDS:
            rows.append((rsid, loci[i], ac[i]))
    return rows


def main() -> int:
    manifest_path = Path(os.environ.get("BV_INPUT_ALIGNED_AC_MANIFEST", "aligned_ac_manifest.txt"))
    union_path = Path(os.environ.get("BV_INPUT_UNION_INDEX", "union_locus_index.json"))

    if not manifest_path.exists():
        print(f"aligned_ac manifest not found: {manifest_path}")
        return 0
    if not union_path.exists():
        print(f"union_index not found: {union_path}")
        return 0

    loci, rsids = load_index(union_path)
    entries = load_manifest(manifest_path)
    if not entries:
        print("no aligned_ac entries found in manifest")
        return 0

    datasets = []
    for ds, path in entries:
        if not path.exists():
            print(f"missing aligned_ac for {ds}: {path}")
            continue
        ac = load_ac(path)
        datasets.append((ds, ac))

    if not datasets:
        print("no aligned_ac files loaded")
        return 0

    print("=== Client aligned_ac summary ===")
    for ds, ac in datasets:
        total = sum(ac)
        nz = sum(1 for x in ac if x)
        mx = max(ac) if ac else 0
        print(f"{ds}: loci={len(ac)} total={total} nonzero={nz} max={mx}")
        print("Top loci:")
        for rank, (val, rsid, locus) in enumerate(top_rows(ac, loci, rsids), start=1):
            print(f"  {rank}\t{val}\t{rsid}\t{locus}")
        print("Overlay rows:")
        for rsid, locus, val in overlay_rows(ac, loci, rsids):
            print(f"  {rsid}\t{locus}\tac={val}")
        print("")

    if len(datasets) >= 2:
        ds1, ac1 = datasets[0]
        ds2, ac2 = datasets[1]
        if len(ac1) == len(ac2):
            diffs = [abs(a - b) for a, b in zip(ac1, ac2)]
            same = all(d == 0 for d in diffs)
            max_diff = max(diffs) if diffs else 0
            mean_diff = sum(diffs) / len(diffs) if diffs else 0
            print("=== Client-to-client diff ===")
            print(f"{ds1} vs {ds2}: identical={same} max_abs_diff={max_diff} mean_abs_diff={mean_diff:.4f}")
        else:
            print(f"aligned_ac length mismatch: {ds1}={len(ac1)} {ds2}={len(ac2)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
