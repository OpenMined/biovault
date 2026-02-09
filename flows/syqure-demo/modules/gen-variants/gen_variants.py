import json
import os
from pathlib import Path

datasite = os.environ.get("BV_CURRENT_DATASITE", "unknown")

counts_by_site = {
    "client1@sandbox.local": {"rs1": 3, "rs2": 1},
    "client2@sandbox.local": {"rs2": 2, "rs3": 4},
    "aggregator@sandbox.local": {"rs1": 1, "rs3": 1, "rs4": 2},
}

counts = counts_by_site.get(datasite, {"rs1": 1, "rs2": 1})

rsids_path = Path(os.environ.get("BV_OUTPUT_RSIDS", "rsids.txt"))
counts_path = Path(os.environ.get("BV_OUTPUT_COUNTS", "counts.json"))

rsids = sorted(counts.keys())
rsids_path.write_text("\n".join(rsids) + "\n", encoding="utf-8")
counts_path.write_text(json.dumps(counts, indent=2, sort_keys=True) + "\n", encoding="utf-8")
