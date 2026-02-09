import json
import os
from pathlib import Path

master_path = Path(os.environ.get("BV_INPUT_MASTER_LIST", "master_list.txt"))
counts_path = Path(os.environ.get("BV_INPUT_COUNTS", "counts.json"))
output_path = Path(os.environ.get("BV_OUTPUT_COUNTS_ARRAY", "counts_array.json"))

master_list = [
    line.strip()
    for line in master_path.read_text(encoding="utf-8").splitlines()
    if line.strip()
]
counts = json.loads(counts_path.read_text(encoding="utf-8"))

aligned = [int(counts.get(rsid, 0)) for rsid in master_list]
output_path.write_text(json.dumps(aligned, indent=2) + "\n", encoding="utf-8")
