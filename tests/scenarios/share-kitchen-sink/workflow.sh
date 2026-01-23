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

PIPELINE_RESULTS_DIR="${BV_RESULTS_DIR}/pipeline"
mkdir -p "$PIPELINE_RESULTS_DIR"

"$BV_BIN" run "$PIPELINE_PATH" \
  --results-dir "$PIPELINE_RESULTS_DIR"
