#!/usr/bin/env bash
set -euo pipefail

INPUT_PATH="${BV_INPUT_COMBINED_PATH:-}"
OUT_PATH="${BV_OUTPUT_COMBINED_COPY:-combined_from_agg.txt}"
TIMEOUT="${BV_INPUT_WAIT_SECONDS:-120}"

if [[ -z "$INPUT_PATH" ]]; then
  echo "Missing combined path: ${INPUT_PATH}" >&2
  exit 1
fi

start="$(date +%s)"
while [[ ! -f "$INPUT_PATH" ]]; do
  if (( $(date +%s) - start > TIMEOUT )); then
    echo "Timed out waiting for ${INPUT_PATH}" >&2
    exit 1
  fi
  sleep 0.2
done

cp -f "$INPUT_PATH" "$OUT_PATH"
