#!/bin/bash
set -euo pipefail

# Read the manifest file and combine all rsids into a sorted unique list
# Manifest format: datasite<TAB>path per line

MANIFEST="${BV_INPUT_RSIDS_MANIFEST:-rsids_manifest.txt}"

# Collect all rsids from manifest paths
rsids=""
while IFS=$'\t' read -r datasite path; do
  [[ -z "$path" ]] && continue
  [[ ! -f "$path" ]] && { echo "Warning: $path not found" >&2; continue; }
  rsids="${rsids}$(cat "$path")"$'\n'
done < "$MANIFEST"

# Sort and deduplicate
echo "$rsids" | sort -u | grep -v '^$' > master_list.txt

echo "Built master list with $(wc -l < master_list.txt | tr -d ' ') unique rsids:"
cat master_list.txt
