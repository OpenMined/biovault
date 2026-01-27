#!/usr/bin/env bash
set -euo pipefail

INPUT_PATH="${BV_INPUT_COMBINED:-}"
OUT_PATH="${BV_OUTPUT_COMBINED:-combined_hello.txt}"

if [[ -z "$INPUT_PATH" || ! -f "$INPUT_PATH" ]]; then
  echo "Missing combined file: ${INPUT_PATH}" >&2
  exit 1
fi

cp -f "$INPUT_PATH" "$OUT_PATH"
