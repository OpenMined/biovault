#!/usr/bin/env bash
set -euo pipefail

UNION_INDEX="${BV_INPUT_UNION_INDEX:-}"
ALIGNED_AC="${BV_INPUT_ALIGNED_AC:-}"
AGG_COUNTS="${BV_INPUT_AGGREGATED_COUNTS:-}"
OUT_JSON="${BV_OUTPUT_REPORT_JSON:-report.json}"
OUT_TSV="${BV_OUTPUT_REPORT_TSV:-report.tsv}"

if [[ -z "$UNION_INDEX" || -z "$ALIGNED_AC" || -z "$AGG_COUNTS" ]]; then
  echo "Missing inputs: BV_INPUT_UNION_INDEX or BV_INPUT_ALIGNED_AC or BV_INPUT_AGGREGATED_COUNTS" >&2
  exit 1
fi

python3 "$(dirname "$0")/report_aggregate.py" \
  --union-index "$UNION_INDEX" \
  --aligned-ac "$ALIGNED_AC" \
  --aggregated "$AGG_COUNTS" \
  --out-json "$OUT_JSON" \
  --out-tsv "$OUT_TSV"

echo "report_json: $OUT_JSON"
echo "report_tsv: $OUT_TSV"
