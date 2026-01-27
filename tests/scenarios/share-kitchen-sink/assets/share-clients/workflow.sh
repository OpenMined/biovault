#!/usr/bin/env bash
set -euo pipefail

INPUT_PATH="${BV_INPUT_HELLO_FILE:-}"
OUT_PATH="${BV_OUTPUT_TAGGED:-tagged_hello.txt}"

if [[ -z "$INPUT_PATH" || ! -f "$INPUT_PATH" ]]; then
  echo "Missing input hello file: ${INPUT_PATH}" >&2
  exit 1
fi

cp -f "$INPUT_PATH" "$OUT_PATH"
printf "tagged_by=%s\n" "${BV_CURRENT_DATASITE:-unknown}" >> "$OUT_PATH"
