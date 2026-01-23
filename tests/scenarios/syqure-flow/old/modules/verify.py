#!/usr/bin/env python3
import json
import os
from pathlib import Path

datasite = os.environ.get("BV_CURRENT_DATASITE", "unknown")
master_path = Path(os.environ.get("BV_INPUT_MASTER_LIST", "master_list.txt"))
my_array_path = Path(os.environ.get("BV_INPUT_MY_ARRAY", "counts_array.json"))
aggregated_path = Path(os.environ.get("BV_INPUT_AGGREGATED", "aggregated_counts.json"))
output_path = Path(os.environ.get("BV_OUTPUT_REPORT", "verification_report.txt"))

master_list = [line.strip() for line in master_path.read_text().splitlines() if line.strip()]
my_array = json.loads(my_array_path.read_text())
aggregated = json.loads(aggregated_path.read_text())

report = []
report.append("=" * 50)
report.append(f"{datasite.upper()} VERIFICATION")
report.append("=" * 50)
report.append(f"Master list:      {master_list}")
report.append(f"My input array:   {my_array}")
report.append(f"Aggregated total: {aggregated}")
report.append(f"My sum:           {sum(my_array)}")
report.append(f"Total sum:        {sum(aggregated)}")
report.append("=" * 50)

report_text = "\n".join(report) + "\n"
output_path.write_text(report_text)
print(report_text)
