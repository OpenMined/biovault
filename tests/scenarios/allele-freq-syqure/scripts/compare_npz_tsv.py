#!/usr/bin/env python3
import argparse
import json
import struct
import zipfile
from pathlib import Path


def load_npz_ac(path: Path):
    try:
        with zipfile.ZipFile(path) as zf:
            if "ac.npy" not in zf.namelist():
                raise RuntimeError(f"ac.npy not found in {path}")
            with zf.open("ac.npy") as f:
                data = f.read()
    except Exception as exc:
        raise RuntimeError(f"failed to read npz: {exc}") from exc

    if not data.startswith(b"\x93NUMPY"):
        prefix = data[:16]
        raise RuntimeError(f"invalid npy header (prefix={prefix!r})")

    major = data[6]
    if major == 1:
        header_len = struct.unpack("<H", data[8:10])[0]
        header_start = 10
    elif major == 2:
        header_len = struct.unpack("<I", data[8:12])[0]
        header_start = 12
    else:
        raise RuntimeError(f"unsupported npy version: {major}")

    header_end = header_start + header_len
    header = data[header_start:header_end].decode("latin1")
    if "fortran_order" in header and "True" in header:
        raise RuntimeError("fortran_order arrays not supported")
    if "'<i8'" not in header and "\"<i8\"" not in header:
        raise RuntimeError("unsupported dtype (expected <i8)")

    shape = ()
    if "shape" in header:
        shape_str = header.split("shape")[1]
        shape_str = shape_str.split(")", 1)[0].split("(", 1)[1]
        shape = tuple(int(x.strip()) for x in shape_str.split(",") if x.strip())
    count = 1
    for dim in shape:
        count *= dim
    fmt = "<" + ("q" * count)
    return list(struct.unpack(fmt, data[header_end:header_end + struct.calcsize(fmt)]))


def load_index(path: Path):
    data = json.loads(path.read_text())
    loci = [str(x) for x in data.get("loci", [])]
    rsids = [str(x) for x in data.get("rsids", [""] * len(loci))]
    if len(rsids) < len(loci):
        rsids.extend([""] * (len(loci) - len(rsids)))
    return loci, rsids


def load_tsv(path: Path):
    rows = {}
    with path.open() as f:
        header = f.readline().strip().split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            locus = parts[idx["locus"]]
            rsid = parts[idx.get("rsid", 1)] if "rsid" in idx else ""
            ac = int(parts[idx["ac"]])
            rows[locus] = (rsid, ac)
    return rows


def check_client(label: str, base: Path):
    npz = base / "allele_freq.npz"
    idx = base / "locus_index.json"
    tsv = base / "allele_freq.tsv"
    if not npz.exists() or not idx.exists() or not tsv.exists():
        print(f"[{label}] missing outputs")
        return None
    try:
        ac = load_npz_ac(npz)
    except Exception as exc:
        print(f"[{label}] npz decode failed: {exc}")
        return None
    loci, rsids = load_index(idx)
    tsv_rows = load_tsv(tsv)
    mismatches = 0
    for i, locus in enumerate(loci[:200]):  # sample 200 to keep it quick
        tsv_ac = tsv_rows.get(locus, ("", None))[1]
        if tsv_ac is None or tsv_ac != ac[i]:
            mismatches += 1
    print(f"[{label}] loci={len(loci)} mismatches(sampled)={mismatches}")
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
