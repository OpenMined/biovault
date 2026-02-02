"""
Build a union locus index from multiple locus_index.json files.
"""

# /// script
# dependencies = ["numpy"]
# ///

import os
from pathlib import Path
import sys

manifest_path = Path(os.environ.get("BV_INPUT_INDEX_MANIFEST", "index_manifest.txt"))
output_path = Path(os.environ.get("BV_OUTPUT_UNION_INDEX", "union_locus_index.json"))
count_path = Path(os.environ.get("BV_OUTPUT_COUNT", "count.txt"))

script_dir = Path(__file__).resolve()
scenarios_dir = script_dir.parents[3]
allele_utils_dir = scenarios_dir / "allele-freq" / "allele-freq"
if allele_utils_dir.exists():
    sys.path.insert(0, str(allele_utils_dir))

from allele_freq_utils import LocusIndex  # type: ignore

indices = []

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
        indices.append(LocusIndex.load(str(file_path)))

if indices:
    union = LocusIndex.from_union(*indices)
else:
    union = LocusIndex.from_union()

union.save(str(output_path))
count_path.write_text(str(len(union)), encoding="utf-8")
