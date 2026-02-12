#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "aggregate: manifest=$BV_INPUT_DATA_MANIFEST"
python3 "$SCRIPT_DIR/aggregate.py" \
  --manifest "$BV_INPUT_DATA_MANIFEST" \
  --out-tsv "$BV_OUTPUT_AGGREGATED_ALLELE_FREQ_TSV" \
  --out-json "$BV_OUTPUT_REPORT_JSON" \
  --out-index "$BV_OUTPUT_UNION_INDEX"
echo "aggregate: done"
