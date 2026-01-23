import json
import os
from pathlib import Path

manifest_path = Path(os.environ.get("BV_INPUT_COUNTS_MANIFEST", "counts_manifest.txt"))
output_path = Path(os.environ.get("BV_OUTPUT_AGGREGATED_COUNTS", "aggregated_counts.json"))

arrays = []
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
    arrays.append(json.loads(file_path.read_text(encoding="utf-8")))

if arrays:
    length = max(len(arr) for arr in arrays)
    totals = [0] * length
    for arr in arrays:
        for idx, value in enumerate(arr):
            totals[idx] += int(value)
else:
    totals = []

output_path.write_text(json.dumps(totals, indent=2) + "\n", encoding="utf-8")
