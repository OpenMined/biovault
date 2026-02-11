"""
Build a union locus index from multiple locus_index.json files.
"""

# /// script
# dependencies = ["numpy"]
# ///

import os
from pathlib import Path
import sys
import numpy as np

manifest_path = Path(os.environ.get("BV_INPUT_INDEX_MANIFEST", "index_manifest.txt"))
output_path = Path(os.environ.get("BV_OUTPUT_UNION_INDEX", "union_locus_index.json"))
count_path = Path(os.environ.get("BV_OUTPUT_COUNT", "count.txt"))
max_loci_raw = os.environ.get("BV_ALLELE_FREQ_MAX_LOCI", "").strip()
force_len_raw = os.environ.get("BV_ALLELE_FREQ_FORCE_ARRAY_LENGTH", "").strip()

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir))

from allele_freq_utils import LocusIndex  # type: ignore

indices = []


def load_index_from_tsv(tsv_path: Path) -> LocusIndex:
    loci = []
    rsids = []
    with tsv_path.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}
        locus_col = col_idx.get("locus_key", col_idx.get("locus"))
        if locus_col is None:
            raise ValueError(f"TSV missing locus/locus_key column: {tsv_path}")
        rsid_col = col_idx.get("rsid")
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if locus_col >= len(parts):
                continue
            loci.append(parts[locus_col])
            rsids.append(parts[rsid_col] if rsid_col is not None and rsid_col < len(parts) else "")
    return LocusIndex(
        loci=np.array(loci, dtype=object),
        rsids=np.array(rsids, dtype=object),
    )

if manifest_path.exists():
    for line in manifest_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        _, path = parts
        path = path.strip()
        if not path:
            continue
        file_path = Path(path)
        if not file_path.exists():
            continue
        if file_path.suffix.lower() == ".json":
            indices.append(LocusIndex.load(str(file_path)))
        else:
            indices.append(load_index_from_tsv(file_path))

if indices:
    union = LocusIndex.from_union(*indices)
else:
    union = LocusIndex.from_union()

if max_loci_raw:
    try:
        max_loci = int(max_loci_raw)
    except ValueError:
        raise SystemExit(f"Invalid BV_ALLELE_FREQ_MAX_LOCI: {max_loci_raw!r}")
    if max_loci > 0 and len(union) > max_loci:
        union = LocusIndex(
            loci=union.loci[:max_loci],
            rsids=union.rsids[:max_loci],
            version=union.version,
        )
        print(f"Capped union loci to {max_loci} via BV_ALLELE_FREQ_MAX_LOCI")

union.save(str(output_path))
if force_len_raw:
    try:
        forced_len = int(force_len_raw)
    except ValueError:
        raise SystemExit(f"Invalid BV_ALLELE_FREQ_FORCE_ARRAY_LENGTH: {force_len_raw!r}")
    if forced_len <= 0:
        raise SystemExit(
            f"BV_ALLELE_FREQ_FORCE_ARRAY_LENGTH must be > 0, got {forced_len}"
        )
    count_path.write_text(str(forced_len), encoding="utf-8")
    print(f"Forced count.txt to {forced_len} via BV_ALLELE_FREQ_FORCE_ARRAY_LENGTH")
else:
    # secure_aggregate consumes concatenated [ac..., an...] vectors from align_counts.
    # The expected MPC array length is therefore 2 * number_of_loci.
    count_path.write_text(str(len(union) * 2), encoding="utf-8")
