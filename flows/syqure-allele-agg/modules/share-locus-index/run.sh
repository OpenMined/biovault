#!/usr/bin/env bash
set -euo pipefail

INPUT="${BV_INPUT_LOCUS_INDEX:-}"
OUTPUT="${BV_OUTPUT_LOCUS_INDEX:-locus_index.tsv}"

if [[ -z "$INPUT" ]]; then
  echo "BV_INPUT_LOCUS_INDEX is required" >&2
  exit 1
fi

cp "$INPUT" "$OUTPUT"
