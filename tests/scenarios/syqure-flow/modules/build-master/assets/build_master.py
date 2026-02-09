import os
from pathlib import Path

manifest_path = Path(os.environ.get("BV_INPUT_RSIDS_MANIFEST", "rsids_manifest.txt"))
output_path = Path(os.environ.get("BV_OUTPUT_MASTER_LIST", "master_list.txt"))

rsids = set()
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
    for rsid in file_path.read_text(encoding="utf-8").splitlines():
        rsid = rsid.strip()
        if rsid:
            rsids.add(rsid)

sorted_rsids = sorted(rsids)
output_path.write_text("\n".join(sorted_rsids) + "\n", encoding="utf-8")

count_path = Path(os.environ.get("BV_OUTPUT_COUNT", "count.txt"))
count_path.write_text(str(len(sorted_rsids)), encoding="utf-8")
