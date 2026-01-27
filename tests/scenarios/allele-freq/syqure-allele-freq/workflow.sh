#!/usr/bin/env bash
set -euo pipefail

BV_BIN="${BV_BIN:-bv}"
PIPELINE_PATH="${BV_ASSETS_DIR}/pipeline.yaml"

if ! command -v "$BV_BIN" >/dev/null 2>&1; then
  echo "bv binary not found: ${BV_BIN}" >&2
  exit 1
fi

if [[ ! -f "$PIPELINE_PATH" ]]; then
  echo "Missing pipeline.yaml at ${PIPELINE_PATH}" >&2
  exit 1
fi

DATA_ROOT="${ALLELE_FREQ_DATA_DIR:-}"
if [[ -z "$DATA_ROOT" ]]; then
  if [[ -n "${BV_SYFTBOX_DATA_DIR:-}" ]]; then
    DATA_ROOT="${BV_SYFTBOX_DATA_DIR}/private/app_data/biovault/allele-freq-data"
  elif [[ -n "${BV_DATASITES_ROOT:-}" ]]; then
    DATA_ROOT="$(dirname "${BV_DATASITES_ROOT}")/private/app_data/biovault/allele-freq-data"
  fi
fi

if [[ -z "$DATA_ROOT" ]]; then
  echo "Unable to resolve allele-freq data root (set ALLELE_FREQ_DATA_DIR)." >&2
  exit 1
fi

SAMPLESHEET="${DATA_ROOT}/samplesheet.csv"
if [[ ! -f "$SAMPLESHEET" ]]; then
  echo "Samplesheet not found: ${SAMPLESHEET}" >&2
  exit 1
fi

PIPELINE_RESULTS_DIR="${BV_RESULTS_DIR}/pipeline"
mkdir -p "$PIPELINE_RESULTS_DIR"

"$BV_BIN" run "$PIPELINE_PATH" \
  --set inputs.samplesheet="$SAMPLESHEET" \
  --results-dir "$PIPELINE_RESULTS_DIR"
