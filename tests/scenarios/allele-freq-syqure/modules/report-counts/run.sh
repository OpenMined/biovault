#!/usr/bin/env bash
set -euo pipefail

UNION_INDEX="${BV_INPUT_UNION_INDEX:-}"
ALIGNED_AC="${BV_INPUT_ALIGNED_AC:-}"

if [[ -z "$UNION_INDEX" || -z "$ALIGNED_AC" ]]; then
  echo "Missing inputs: BV_INPUT_UNION_INDEX or BV_INPUT_ALIGNED_AC" >&2
  exit 1
fi

python3 "$(dirname "$0")/report_counts.py" \
  --union-index "$UNION_INDEX" \
  --aligned-ac "$ALIGNED_AC" \
  ${BV_OUTPUT_REPORT_JSON:+--out-json "$BV_OUTPUT_REPORT_JSON"} \
  ${BV_OUTPUT_REPORT_TSV:+--out-tsv "$BV_OUTPUT_REPORT_TSV"}
