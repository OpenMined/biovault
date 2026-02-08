#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


APOL1 = {"rs60910145", "rs71785313", "rs73885319"}
THAL = {"rs33985472", "rs33971634"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--union-index", required=True)
    parser.add_argument("--aligned-ac", required=True)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--out-json")
    parser.add_argument("--out-tsv")
    args = parser.parse_args()

    union_index = Path(args.union_index)
    aligned_ac = Path(args.aligned_ac)
    if not union_index.exists() or not aligned_ac.exists():
        raise SystemExit("missing union_index or aligned_ac")

    index = json.loads(union_index.read_text())
    loci = index.get("loci", [])
    rsids = index.get("rsids", [""] * len(loci))
    ac = json.loads(aligned_ac.read_text())

    n = len(ac)
    total = sum(ac)
    nonzero = sum(1 for x in ac if x)
    max_val = max(ac) if ac else 0

    print(f"aligned_ac: {aligned_ac}")
    print(f"loci={len(loci)} ac_len={n} total={total} nonzero={nonzero} max={max_val}")

    top_n = min(args.top, n)
    top_idx = sorted(range(n), key=lambda i: ac[i], reverse=True)[:top_n]
    print("\nTop loci by count:")
    print("rank\tac\trsid\tlocus")
    for rank, i in enumerate(top_idx, start=1):
        if ac[i] == 0:
            continue
        rs = rsids[i] if i < len(rsids) else ""
        locus = loci[i] if i < len(loci) else ""
        print(f"{rank}\t{ac[i]}\t{rs}\t{locus}")

    def report(name, rs_set):
        rows = []
        for locus, rsid, count in zip(loci, rsids, ac):
            if rsid in rs_set:
                rows.append((rsid, locus, count))
        print(f"\n{name} rows:")
        if not rows:
            print("  (none)")
            return
        for rsid, locus, count in rows:
            print(f"  {rsid}\t{locus}\tac={count}")

    report("APOL1", APOL1)
    report("THAL", THAL)

    if args.out_tsv:
        out_path = Path(args.out_tsv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8") as f:
            f.write("rank\tac\trsid\tlocus\n")
            for rank, i in enumerate(top_idx, start=1):
                if ac[i] == 0:
                    continue
                rs = rsids[i] if i < len(rsids) else ""
                locus = loci[i] if i < len(loci) else ""
                f.write(f"{rank}\t{ac[i]}\t{rs}\t{locus}\n")

    if args.out_json:
        out_path = Path(args.out_json)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        def collect(rs_set):
            rows = []
            for locus, rsid, count in zip(loci, rsids, ac):
                if rsid in rs_set:
                    rows.append({"rsid": rsid, "locus": locus, "ac": int(count)})
            return rows
        top_rows = []
        for rank, i in enumerate(top_idx, start=1):
            if ac[i] == 0:
                continue
            rs = rsids[i] if i < len(rsids) else ""
            locus = loci[i] if i < len(loci) else ""
            top_rows.append({"rank": rank, "ac": int(ac[i]), "rsid": rs, "locus": locus})
        payload = {
            "aligned_ac": str(aligned_ac),
            "loci": len(loci),
            "ac_len": n,
            "total": int(total),
            "nonzero": int(nonzero),
            "max": int(max_val),
            "top": top_rows,
            "apol1": collect(APOL1),
            "thal": collect(THAL),
        }
        out_path.write_text(json.dumps(payload, separators=(",", ":")) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
