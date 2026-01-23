#!/usr/bin/env bash
set -euo pipefail

OUT_PATH="${BV_OUTPUT_HELLO:-hello.txt}"

printf "current=%s\n" "${BV_CURRENT_DATASITE:-unknown}" > "$OUT_PATH"
printf "index=%s\n" "${BV_DATASITE_INDEX:-unknown}" >> "$OUT_PATH"
printf "datasites=%s\n" "${BV_DATASITES:-}" >> "$OUT_PATH"
printf "run_id=%s\n" "${BV_RUN_ID:-}" >> "$OUT_PATH"
