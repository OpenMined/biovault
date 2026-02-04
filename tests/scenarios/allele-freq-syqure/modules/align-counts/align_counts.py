"""
Align allele frequency arrays to the union index and emit combined AC+AN vector.

The output is a concatenated vector: [ac_0, ac_1, ..., ac_n, an_0, an_1, ..., an_n]
This allows a single MPC aggregation instead of two separate ones.
"""

# /// script
# dependencies = ["numpy"]
# ///

import json
import os
from pathlib import Path
import sys

tsv_path = os.environ.get("BV_INPUT_ALLELE_FREQ_TSV", "")
union_path = os.environ.get("BV_INPUT_UNION_INDEX", "")

out_counts = Path(os.environ.get("BV_OUTPUT_ALIGNED_COUNTS", "aligned_counts.json"))

if not union_path or not Path(union_path).exists():
    raise SystemExit("Missing union index file")

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir))

from allele_freq_utils import AlleleFreqData, LocusIndex  # type: ignore

union = LocusIndex.load(union_path)
n = len(union)

if tsv_path and Path(tsv_path).exists():
    data = AlleleFreqData.from_tsv(tsv_path)
    aligned = data.align_to(union)
    ac = aligned.ac.tolist()
    an = aligned.an.tolist()
else:
    ac = [0 for _ in range(n)]
    an = [0 for _ in range(n)]

# Concatenate AC and AN into single vector: [ac..., an...]
combined = ac + an

out_counts.write_text(json.dumps(combined, separators=(",", ":")) + "\n", encoding="utf-8")
print(f"Aligned counts: {n} loci, {len(combined)} total elements (AC+AN concatenated)")
