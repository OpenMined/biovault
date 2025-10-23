#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BV_BIN="${ROOT_DIR}/bv"

DATA_DIR="${ROOT_DIR}/cli/examples/pipeline/data"
FILTER_PROJECT="${ROOT_DIR}/cli/examples/pipeline/filter-samples"
COUNT_PROJECT="${ROOT_DIR}/cli/examples/pipeline/count-lines"
SUM_PROJECT="${ROOT_DIR}/cli/examples/pipeline/sum-total"

# Stage 1: filter samplesheet
"${BV_BIN}" run "${FILTER_PROJECT}" \
  --results-dir results_filter \
  --samplesheet "${DATA_DIR}/participants.csv" \
  --data_dir "${DATA_DIR}" \
  --param.extension .txt

FILTER_RESULTS_DIR="${FILTER_PROJECT}/results_filter"

# Stage 2: annotate line counts
"${BV_BIN}" run "${COUNT_PROJECT}" \
  --results-dir results_counts \
  --samplesheet "${FILTER_RESULTS_DIR}/filtered_samplesheet.csv" \
  --data_dir "${DATA_DIR}"

COUNT_RESULTS_DIR="${COUNT_PROJECT}/results_counts"

# Stage 3: sum total line counts
"${BV_BIN}" run "${SUM_PROJECT}" \
  --results-dir results_sum \
  --samplesheet "${COUNT_RESULTS_DIR}/line_counts.csv" \
  --data_dir "${DATA_DIR}"
