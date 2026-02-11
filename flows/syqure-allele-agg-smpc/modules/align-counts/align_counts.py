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
force_len_raw = os.environ.get("BV_ALLELE_FREQ_FORCE_ARRAY_LENGTH", "").strip()

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

if force_len_raw:
    try:
        force_len = int(force_len_raw)
    except ValueError:
        raise SystemExit(f"Invalid BV_ALLELE_FREQ_FORCE_ARRAY_LENGTH: {force_len_raw!r}")
    if force_len <= 0:
        raise SystemExit(
            f"BV_ALLELE_FREQ_FORCE_ARRAY_LENGTH must be > 0, got {force_len}"
        )
    # Prefer non-zero values for debug visibility so tiny-vector smoke tests
    # still exercise arithmetic instead of all-zero paths.
    non_zero = [v for v in combined if v != 0]
    if len(non_zero) >= force_len:
        combined = non_zero[:force_len]
    else:
        combined = non_zero + [0 for _ in range(force_len - len(non_zero))]
    print(
        f"Forced aligned counts vector length to {force_len} "
        f"via BV_ALLELE_FREQ_FORCE_ARRAY_LENGTH"
    )

out_counts.write_text(json.dumps(combined, separators=(",", ":")) + "\n", encoding="utf-8")
print(f"Aligned counts: {n} loci, {len(combined)} total elements (AC+AN concatenated)")

# Compute and output AN bounds for Chebyshev interval
an_positive = [v for v in an if v > 0]
an_min = min(an_positive) if an_positive else 0
an_max = max(an_positive) if an_positive else 0
out_an_bounds = Path(os.environ.get("BV_OUTPUT_AN_BOUNDS", "an_bounds.txt"))
out_an_bounds.write_text(f"{an_min}\n{an_max}\n", encoding="utf-8")
print(f"AN bounds: min={an_min}, max={an_max} (from {len(an_positive)} positive values)")
