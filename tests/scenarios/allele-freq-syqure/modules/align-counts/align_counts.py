"""
Align allele frequency arrays to the union index and emit AC array.
"""

# /// script
# dependencies = ["numpy"]
# ///

import json
import os
from pathlib import Path
import sys

npz_path = os.environ.get("BV_INPUT_ALLELE_FREQ_NPZ", "")
index_path = os.environ.get("BV_INPUT_LOCUS_INDEX", "")
union_path = os.environ.get("BV_INPUT_UNION_INDEX", "")

out_ac = Path(os.environ.get("BV_OUTPUT_ALIGNED_AC", "aligned_ac.json"))

if not union_path or not Path(union_path).exists():
    raise SystemExit("Missing union index file")

script_dir = Path(__file__).resolve()
scenarios_dir = script_dir.parents[3]
allele_utils_dir = scenarios_dir / "allele-freq" / "allele-freq"
if allele_utils_dir.exists():
    sys.path.insert(0, str(allele_utils_dir))

from allele_freq_utils import AlleleFreqData, LocusIndex  # type: ignore

union = LocusIndex.load(union_path)

if npz_path and index_path and Path(npz_path).exists() and Path(index_path).exists():
    data = AlleleFreqData.load(npz_path, index_path)
    aligned = data.align_to(union)
    ac = aligned.ac.tolist()
else:
    ac = [0 for _ in range(len(union))]

out_ac.write_text(json.dumps(ac, separators=(",", ":")) + "\n", encoding="utf-8")
