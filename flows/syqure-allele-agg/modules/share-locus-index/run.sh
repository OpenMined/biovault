#!/usr/bin/env bash
set -euo pipefail

INPUT="${BV_INPUT_LOCUS_INDEX:-}"
OUTPUT="${BV_OUTPUT_LOCUS_INDEX:-locus_index.tsv}"

if [[ -z "$INPUT" ]]; then
  echo "BV_INPUT_LOCUS_INDEX is required" >&2
  exit 1
fi

# Share only the locus key column (column 1) for build_master.
awk -F'\t' 'NF { sub(/\r$/, "", $1); print $1 }' "$INPUT" > "$OUTPUT"
