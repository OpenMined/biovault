#!/usr/bin/env bash
set -euo pipefail

MANIFEST_PATH="${BV_INPUT_SHARED_PATHS:-}"
OUT_PATH="${BV_OUTPUT_COMBINED:-combined_hello.txt}"
TIMEOUT="${BV_INPUT_WAIT_SECONDS:-120}"

if [[ -z "$MANIFEST_PATH" || ! -f "$MANIFEST_PATH" ]]; then
  echo "Missing shared_paths manifest: ${MANIFEST_PATH}" >&2
  exit 1
fi

: > "$OUT_PATH"

wait_for_file() {
  local path="$1"
  local start
  start="$(date +%s)"
  while [[ ! -f "$path" ]]; do
    if (( $(date +%s) - start > TIMEOUT )); then
      echo "Timed out waiting for ${path}" >&2
      return 1
    fi
    sleep 0.2
  done
}

while IFS=$'\t' read -r site path; do
  if [[ -z "$site" || -z "$path" ]]; then
    continue
  fi
  wait_for_file "$path"
  content=$(tr '\n' '|' < "$path" | sed 's/|$//')
  printf "%s\t%s\n" "$site" "$content" >> "$OUT_PATH"
done < "$MANIFEST_PATH"
