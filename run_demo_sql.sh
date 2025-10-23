#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_FILE="${ROOT_DIR}/pipeline_sql.yaml"
RESULTS_DIR="${ROOT_DIR}/results/demo-sql"
SAMPLESHEET="${ROOT_DIR}/cli/examples/pipeline/data/participants.csv"
DATA_DIR="${ROOT_DIR}/cli/examples/pipeline/data"
BV="${ROOT_DIR}/bv"

"${BV}" pipeline validate --diagram "${PIPELINE_FILE}"

"${BV}" run "${PIPELINE_FILE}" \
  --set inputs.samplesheet="${SAMPLESHEET}" \
  --set inputs.data_dir="${DATA_DIR}" \
  --results-dir "${RESULTS_DIR}"
