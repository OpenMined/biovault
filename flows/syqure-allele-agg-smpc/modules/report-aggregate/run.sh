#!/usr/bin/env bash
set -euo pipefail

UNION_INDEX="${BV_INPUT_UNION_INDEX:-}"
LOCAL_COUNTS="${BV_INPUT_LOCAL_COUNTS:-}"
AGG_COUNTS="${BV_INPUT_AGGREGATED_COUNTS:-}"
OUT_JSON="${BV_OUTPUT_REPORT_JSON:-report.json}"
OUT_TSV="${BV_OUTPUT_REPORT_TSV:-report.tsv}"
OUT_AGG_TSV="${BV_OUTPUT_AGGREGATED_ALLELE_FREQ_TSV:-aggregated_allele_freq.tsv}"

if [[ -z "$UNION_INDEX" || -z "$LOCAL_COUNTS" || -z "$AGG_COUNTS" ]]; then
  echo "Missing inputs:" >&2
  echo "  BV_INPUT_UNION_INDEX=$UNION_INDEX" >&2
  echo "  BV_INPUT_LOCAL_COUNTS=$LOCAL_COUNTS" >&2
  echo "  BV_INPUT_AGGREGATED_COUNTS=$AGG_COUNTS" >&2
  exit 1
fi

python3 "$(dirname "$0")/report_aggregate.py" \
  --union-index "$UNION_INDEX" \
  --local-counts "$LOCAL_COUNTS" \
  --aggregated-counts "$AGG_COUNTS" \
  --out-json "$OUT_JSON" \
  --out-tsv "$OUT_TSV" \
  --out-agg-tsv "$OUT_AGG_TSV"
