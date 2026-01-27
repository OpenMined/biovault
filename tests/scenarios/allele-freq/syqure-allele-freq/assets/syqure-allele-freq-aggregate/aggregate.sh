#!/usr/bin/env bash
set -euo pipefail

PID="${BV_DATASITE_INDEX:-}"
CURRENT_EMAIL="${BV_CURRENT_DATASITE:-}"
PARTY_EMAILS="${BV_DATASITES:-}"
RUN_ID="${BV_RUN_ID:-}"
PROJECT_DIR="${BV_PROJECT_DIR:-}"
ASSETS_DIR="${BV_ASSETS_DIR:-}"
DATASITES_ROOT="${BV_DATASITES_ROOT:-}"
RESULTS_DIR="${BV_RESULTS_DIR:-}"
POLL_MS="${SEQURE_FILE_POLL_MS:-50}"
WAIT_TIMEOUT="${SEQURE_SYFT_WAIT_TIMEOUT:-300}"
WAIT_STRICT="${SEQURE_SYFT_WAIT_STRICT:-0}"

MODE="${SEQURE_MODE:-native}"
NATIVE_BIN="${SEQURE_NATIVE_BIN:-}"

INPUT_NPZ="${BV_INPUT_ALLELE_FREQ_NPZ:-}"
INPUT_INDEX="${BV_INPUT_LOCUS_INDEX:-}"

if [[ -z "$PID" || -z "$CURRENT_EMAIL" || -z "$PARTY_EMAILS" ]]; then
  echo "Missing required datasite context (BV_DATASITE_INDEX, BV_CURRENT_DATASITE, BV_DATASITES)" >&2
  exit 1
fi

if [[ -z "$RUN_ID" ]]; then
  RUN_ID="$(date +%Y%m%d%H%M%S)"
fi

if [[ -z "$DATASITES_ROOT" ]]; then
  echo "Missing BV_DATASITES_ROOT" >&2
  exit 1
fi

if [[ -z "$RESULTS_DIR" ]]; then
  RESULTS_DIR="${PROJECT_DIR}/results"
fi

mkdir -p "$RESULTS_DIR"

if [[ "$MODE" != "native" ]]; then
  echo "Only SEQURE_MODE=native is supported for syqure-allele-freq" >&2
  exit 1
fi

if [[ -z "$NATIVE_BIN" ]]; then
  for candidate in \
    "${PROJECT_DIR}/../../../syqure/target/release/syqure" \
    "${PROJECT_DIR}/../../../syqure/target/debug/syqure" \
    "$(dirname "$(dirname "${BV_SYFTBOX_DATA_DIR:-/tmp}")")/syqure/target/release/syqure" \
    "$(dirname "$(dirname "${BV_SYFTBOX_DATA_DIR:-/tmp}")")/syqure/target/debug/syqure" \
    "$(dirname "$(dirname "$(dirname "${BV_SYFTBOX_DATA_DIR:-/tmp}")")")/syqure/target/release/syqure" \
    "$(dirname "$(dirname "$(dirname "${BV_SYFTBOX_DATA_DIR:-/tmp}")")")/syqure/target/debug/syqure"; do
    if [[ -x "$candidate" ]]; then
      NATIVE_BIN="$candidate"
      break
    fi
  done
fi

if [[ -z "$NATIVE_BIN" || ! -x "$NATIVE_BIN" ]]; then
  echo "Native syqure binary not found. Set SEQURE_NATIVE_BIN." >&2
  exit 1
fi

echo "[PID=$PID] Using native syqure: $NATIVE_BIN" >&2

IFS=',' read -r -a PARTY_LIST <<< "$PARTY_EMAILS"
CP_COUNT="${#PARTY_LIST[@]}"

if [[ "$PID" -lt 0 || "$PID" -ge "$CP_COUNT" ]]; then
  echo "PID ${PID} out of range for party list (${CP_COUNT})" >&2
  exit 1
fi

LOCAL_BASE="${DATASITES_ROOT}/${CURRENT_EMAIL}/shared/syqure/${RUN_ID}"

write_acl() {
  local dir="$1"
  local peer_a="$2"
  local peer_b="$3"
  mkdir -p "$dir"
  cat > "${dir}/syft.pub.yaml" <<EOF_ACL
rules:
  - pattern: "**"
    access:
      admin: []
      read:
        - ${peer_a}
        - ${peer_b}
      write:
        - ${peer_a}
        - ${peer_b}
EOF_ACL
}

wait_for() {
  local path="$1"
  local timeout="${2:-120}"
  local start
  start="$(date +%s)"
  while [[ ! -f "$path" ]]; do
    if (( $(date +%s) - start > timeout )); then
      echo "Timed out waiting for ${path}" >&2
      if [[ "$WAIT_STRICT" == "1" ]]; then
        exit 1
      fi
      return 0
    fi
    sleep 0.2
  done
}

for ((i=0; i<CP_COUNT; i++)); do
  if [[ "$i" -eq "$PID" ]]; then
    continue
  fi
  peer_email="${PARTY_LIST[$i]}"
  write_acl "${LOCAL_BASE}/${PID}_to_${i}" "${CURRENT_EMAIL}" "${peer_email}"
  wait_for "${DATASITES_ROOT}/${peer_email}/shared/syqure/${RUN_ID}/${i}_to_${PID}/syft.pub.yaml" "${WAIT_TIMEOUT}"
done

AGG_EMAIL="${PARTY_LIST[0]}"
LOCAL_INDEX_PATH="${RESULTS_DIR}/locus_index.json"

if [[ -n "$INPUT_INDEX" && -f "$INPUT_INDEX" ]]; then
  cp -f "$INPUT_INDEX" "$LOCAL_INDEX_PATH"
else
  cat > "$LOCAL_INDEX_PATH" <<'EOF_INDEX'
{"version":"1.0","n_loci":0,"loci":[],"rsids":[]}
EOF_INDEX
fi

if [[ "$PID" -ne 0 ]]; then
  cp -f "$LOCAL_INDEX_PATH" "${LOCAL_BASE}/${PID}_to_0/locus_index.json"
fi

UNION_INDEX_PATH="${RESULTS_DIR}/union_locus_index.json"

if [[ "$PID" -eq 0 ]]; then
  for ((i=1; i<CP_COUNT; i++)); do
    peer_email="${PARTY_LIST[$i]}"
    wait_for "${DATASITES_ROOT}/${peer_email}/shared/syqure/${RUN_ID}/${i}_to_0/locus_index.json" "${WAIT_TIMEOUT}"
  done

  UNION_INDEX_OUT="${UNION_INDEX_PATH}" PYTHONPATH="${ASSETS_DIR}" python3 - \
    "${LOCAL_INDEX_PATH}" \
    "${DATASITES_ROOT}/${PARTY_LIST[1]}/shared/syqure/${RUN_ID}/1_to_0/locus_index.json" \
    "${DATASITES_ROOT}/${PARTY_LIST[2]}/shared/syqure/${RUN_ID}/2_to_0/locus_index.json" <<'PY'
import os
import sys
import numpy as np
from allele_freq_utils import LocusIndex

index_paths = sys.argv[1:]
indices = []
for path in index_paths:
    if os.path.isfile(path):
        indices.append(LocusIndex.load(path))

if indices:
    union = LocusIndex.from_union(*indices)
else:
    union = LocusIndex(loci=np.array([], dtype=object), rsids=np.array([], dtype=object))

union.save(os.environ["UNION_INDEX_OUT"])
PY

  cp -f "$UNION_INDEX_PATH" "${LOCAL_BASE}/0_to_1/union_locus_index.json"
  cp -f "$UNION_INDEX_PATH" "${LOCAL_BASE}/0_to_2/union_locus_index.json"
else
  wait_for "${DATASITES_ROOT}/${AGG_EMAIL}/shared/syqure/${RUN_ID}/0_to_${PID}/union_locus_index.json" "${WAIT_TIMEOUT}"
  cp -f "${DATASITES_ROOT}/${AGG_EMAIL}/shared/syqure/${RUN_ID}/0_to_${PID}/union_locus_index.json" "$UNION_INDEX_PATH"
fi

ALIGNED_DIR="${RESULTS_DIR}/aligned"
mkdir -p "$ALIGNED_DIR"
ALIGNED_AC="${ALIGNED_DIR}/aligned_ac.txt"
ALIGNED_AN="${ALIGNED_DIR}/aligned_an.txt"

PYTHONPATH="${ASSETS_DIR}" python3 - <<'PY'
import os
import numpy as np
from allele_freq_utils import AlleleFreqData, LocusIndex

npz_path = os.environ.get("INPUT_NPZ", "")
index_path = os.environ.get("INPUT_INDEX", "")
union_path = os.environ["UNION_INDEX"]
ac_out = os.environ["ALIGNED_AC"]
an_out = os.environ["ALIGNED_AN"]

union = LocusIndex.load(union_path)

if npz_path and index_path and os.path.isfile(npz_path) and os.path.isfile(index_path):
    data = AlleleFreqData.load(npz_path, index_path)
    aligned = data.align_to(union)
    ac = aligned.ac
    an = aligned.an
else:
    n = len(union.loci)
    ac = np.zeros(n, dtype=np.int64)
    an = np.zeros(n, dtype=np.int64)

np.savetxt(ac_out, ac, fmt="%d")
np.savetxt(an_out, an, fmt="%d")
PY \
  INPUT_NPZ="${INPUT_NPZ}" \
  INPUT_INDEX="${INPUT_INDEX}" \
  UNION_INDEX="${UNION_INDEX_PATH}" \
  ALIGNED_AC="${ALIGNED_AC}" \
  ALIGNED_AN="${ALIGNED_AN}"

PROGRAM_PATH="${ASSETS_DIR}/allele_freq_sum.codon"
if [[ ! -f "$PROGRAM_PATH" ]]; then
  echo "Missing codon program: ${PROGRAM_PATH}" >&2
  exit 1
fi

export SEQURE_TRANSPORT=file
export SEQURE_FILE_DIR="shared/syqure/${RUN_ID}"
export SEQURE_FILE_POLL_MS="${POLL_MS}"
export SEQURE_CP_COUNT="${CP_COUNT}"
export SEQURE_PARTY_EMAILS="${PARTY_EMAILS}"
export SEQURE_DATASITES_ROOT="${DATASITES_ROOT}"
export SEQURE_LOCAL_EMAIL="${CURRENT_EMAIL}"

export BV_INPUT_AC="${ALIGNED_AC}"
export BV_INPUT_AN="${ALIGNED_AN}"
export BV_OUTPUT_SUM_AC="${RESULTS_DIR}/sum_ac.txt"
export BV_OUTPUT_SUM_AN="${RESULTS_DIR}/sum_an.txt"

LOG_DIR="${SEQURE_LOG_DIR:-}"
LOG_PATH=""
if [[ -n "$LOG_DIR" ]]; then
  mkdir -p "${RESULTS_DIR}/${LOG_DIR}"
  LOG_PATH="${RESULTS_DIR}/${LOG_DIR}/syqure-${RUN_ID}-pid${PID}.log"
fi

if [[ -n "$LOG_PATH" ]]; then
  "$NATIVE_BIN" "$PROGRAM_PATH" -- "$PID" 2>&1 | tee "$LOG_PATH"
  exit "${PIPESTATUS[0]}"
fi

exec "$NATIVE_BIN" "$PROGRAM_PATH" -- "$PID"
