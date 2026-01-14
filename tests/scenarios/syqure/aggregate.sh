#!/usr/bin/env bash
set -euo pipefail

PID="${BV_DATASITE_INDEX:-}"
CURRENT_EMAIL="${BV_CURRENT_DATASITE:-}"
PARTY_EMAILS="${BV_DATASITES:-}"
RUN_ID="${BV_RUN_ID:-}"
PROJECT_DIR="${BV_PROJECT_DIR:-}"
ASSETS_DIR="${BV_ASSETS_DIR:-}"
DATASITES_ROOT="${BV_DATASITES_ROOT:-}"
POLL_MS="${SEQURE_FILE_POLL_MS:-50}"
WAIT_TIMEOUT="${SEQURE_SYFT_WAIT_TIMEOUT:-300}"
WAIT_STRICT="${SEQURE_SYFT_WAIT_STRICT:-0}"

# Mode: native (default) or docker
MODE="${SEQURE_MODE:-native}"
NATIVE_BIN="${SEQURE_NATIVE_BIN:-}"
IMAGE_NAME="${SEQURE_IMAGE_NAME:-syqure-cli-files}"
PLATFORM="${SEQURE_PLATFORM:-linux/amd64}"

if [[ -z "$PID" || -z "$CURRENT_EMAIL" || -z "$PARTY_EMAILS" ]]; then
  echo "Missing required datasite context (BV_DATASITE_INDEX, BV_CURRENT_DATASITE, BV_DATASITES)" >&2
  exit 1
fi

if [[ -z "$RUN_ID" ]]; then
  RUN_ID="$(date +%Y%m%d%H%M%S)"
fi

if [[ -z "$DATASITES_ROOT" ]]; then
  if [[ -n "${SYFTBOX_DATA_DIR:-}" ]]; then
    DATASITES_ROOT="${SYFTBOX_DATA_DIR}/datasites"
  fi
fi

if [[ -z "$DATASITES_ROOT" || -z "$PROJECT_DIR" ]]; then
  echo "Missing required paths (BV_PROJECT_DIR or BV_DATASITES_ROOT)" >&2
  exit 1
fi

mkdir -p "$DATASITES_ROOT"

# Find native binary if in native mode
if [[ "$MODE" == "native" ]]; then
  if [[ -z "$NATIVE_BIN" ]]; then
    # Try common locations relative to project
    # SYFTBOX_DATA_DIR is typically sandbox/{client}, so we need ../../syqure (workspace2/syqure)
    for candidate in \
      "${PROJECT_DIR}/../../../syqure/target/release/syqure" \
      "${PROJECT_DIR}/../../../syqure/target/debug/syqure" \
      "$(dirname "$(dirname "${SYFTBOX_DATA_DIR:-/tmp}")")/syqure/target/release/syqure" \
      "$(dirname "$(dirname "${SYFTBOX_DATA_DIR:-/tmp}")")/syqure/target/debug/syqure" \
      "$(dirname "$(dirname "$(dirname "${SYFTBOX_DATA_DIR:-/tmp}")")")/syqure/target/release/syqure" \
      "$(dirname "$(dirname "$(dirname "${SYFTBOX_DATA_DIR:-/tmp}")")")/syqure/target/debug/syqure"; do
      if [[ -x "$candidate" ]]; then
        NATIVE_BIN="$candidate"
        break
      fi
    done
  fi
  if [[ -z "$NATIVE_BIN" || ! -x "$NATIVE_BIN" ]]; then
    echo "Native syqure binary not found. Set SEQURE_NATIVE_BIN or use SEQURE_MODE=docker" >&2
    exit 1
  fi
  echo "[PID=$PID] Using native syqure: $NATIVE_BIN" >&2
fi

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
  cat > "${dir}/syft.pub.yaml" <<EOF
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
EOF
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

# Create ACL directories for outgoing channels
for ((i=0; i<CP_COUNT; i++)); do
  if [[ "$i" -eq "$PID" ]]; then
    continue
  fi
  peer_email="${PARTY_LIST[$i]}"
  write_acl "${LOCAL_BASE}/${PID}_to_${i}" "${CURRENT_EMAIL}" "${peer_email}"
done

# Wait for peer ACLs to sync
for ((i=0; i<CP_COUNT; i++)); do
  if [[ "$i" -eq "$PID" ]]; then
    continue
  fi
  peer_email="${PARTY_LIST[$i]}"
  wait_for "${DATASITES_ROOT}/${peer_email}/shared/syqure/${RUN_ID}/${i}_to_${PID}/syft.pub.yaml" "${WAIT_TIMEOUT}"
done

# Copy input file if needed
if [[ -n "$ASSETS_DIR" && -f "${ASSETS_DIR}/${PID}.txt" && ! -f "${PROJECT_DIR}/${PID}.txt" ]]; then
  cp -f "${ASSETS_DIR}/${PID}.txt" "${PROJECT_DIR}/${PID}.txt"
fi

# Resolve input path
HOST_INPUT_PATH="${BV_INPUT_MY_VAL:-}"
if [[ -z "$HOST_INPUT_PATH" ]]; then
  if [[ -f "${PROJECT_DIR}/${PID}.txt" ]]; then
    HOST_INPUT_PATH="${PROJECT_DIR}/${PID}.txt"
  elif [[ -f "${ASSETS_DIR}/${PID}.txt" ]]; then
    HOST_INPUT_PATH="${ASSETS_DIR}/${PID}.txt"
  fi
fi

# Resolve output path
HOST_OUTPUT_PATH=""
if [[ -n "${BV_OUTPUT_AGGREGATE_RESULT:-}" ]]; then
  if [[ "${BV_OUTPUT_AGGREGATE_RESULT}" == "$PROJECT_DIR/"* ]]; then
    HOST_OUTPUT_PATH="${BV_OUTPUT_AGGREGATE_RESULT}"
  elif [[ "${BV_OUTPUT_AGGREGATE_RESULT}" != /* ]]; then
    HOST_OUTPUT_PATH="${PROJECT_DIR}/${BV_OUTPUT_AGGREGATE_RESULT}"
  fi
fi

if [[ -z "$HOST_OUTPUT_PATH" && -n "${BV_RESULTS_DIR:-}" ]]; then
  if [[ "${BV_RESULTS_DIR}" == "$PROJECT_DIR/"* ]]; then
    HOST_OUTPUT_PATH="${BV_RESULTS_DIR}/${PID}_out.txt"
  elif [[ "${BV_RESULTS_DIR}" != /* ]]; then
    HOST_OUTPUT_PATH="${PROJECT_DIR}/${BV_RESULTS_DIR}/${PID}_out.txt"
  fi
fi

if [[ -n "$HOST_OUTPUT_PATH" ]]; then
  mkdir -p "$(dirname "${HOST_OUTPUT_PATH}")"
fi

# Find program path
PROGRAM_PATH="${PROJECT_DIR}/assets/aggregate.codon"
if [[ ! -f "$PROGRAM_PATH" && -f "${PROJECT_DIR}/aggregate.codon" ]]; then
  PROGRAM_PATH="${PROJECT_DIR}/aggregate.codon"
fi

LOG_DIR="${SEQURE_LOG_DIR:-}"
LOG_PATH=""
if [[ -n "$LOG_DIR" ]]; then
  mkdir -p "${PROJECT_DIR}/${LOG_DIR}"
  LOG_PATH="${PROJECT_DIR}/${LOG_DIR}/syqure-${RUN_ID}-pid${PID}.log"
fi

# Run syqure (native or docker)
if [[ "$MODE" == "native" ]]; then
  # Native mode
  export SEQURE_TRANSPORT=file
  export SEQURE_FILE_DIR="shared/syqure/${RUN_ID}"
  export SEQURE_FILE_POLL_MS="${POLL_MS}"
  export SEQURE_CP_COUNT="${CP_COUNT}"
  export SEQURE_PARTY_EMAILS="${PARTY_EMAILS}"
  export SEQURE_DATASITES_ROOT="${DATASITES_ROOT}"
  export SEQURE_LOCAL_EMAIL="${CURRENT_EMAIL}"
  ${SEQURE_FILE_KEEP:+export SEQURE_FILE_KEEP="${SEQURE_FILE_KEEP}"}
  ${SEQURE_FILE_DEBUG:+export SEQURE_FILE_DEBUG="${SEQURE_FILE_DEBUG}"}
  ${HOST_INPUT_PATH:+export BV_INPUT_MY_VAL="${HOST_INPUT_PATH}"}
  ${HOST_OUTPUT_PATH:+export BV_OUTPUT_AGGREGATE_RESULT="${HOST_OUTPUT_PATH}"}

  if [[ -n "$LOG_PATH" ]]; then
    "$NATIVE_BIN" "$PROGRAM_PATH" -- "$PID" 2>&1 | tee "$LOG_PATH"
    exit "${PIPESTATUS[0]}"
  else
    exec "$NATIVE_BIN" "$PROGRAM_PATH" -- "$PID"
  fi
else
  # Docker mode
  DATASITES_ROOT_IN_CONTAINER="/datasites"
  MOUNT_SPEC="-v ${DATASITES_ROOT}:${DATASITES_ROOT_IN_CONTAINER}"
  if ! docker run --rm ${MOUNT_SPEC} alpine:3.19 sh -c 'test -d /datasites' >/dev/null 2>&1; then
    parent_dir="$(dirname "$DATASITES_ROOT")"
    DATASITES_ROOT_IN_CONTAINER="/datasites-root/datasites"
    MOUNT_SPEC="-v ${parent_dir}:/datasites-root"
  fi

  INPUT_ENV=""
  if [[ -n "$HOST_INPUT_PATH" && "$HOST_INPUT_PATH" == "$PROJECT_DIR/"* ]]; then
    rel="${HOST_INPUT_PATH#$PROJECT_DIR/}"
    INPUT_ENV="/workspace/project/${rel}"
  fi

  OUTPUT_ENV=""
  if [[ -n "$HOST_OUTPUT_PATH" && "$HOST_OUTPUT_PATH" == "$PROJECT_DIR/"* ]]; then
    rel="${HOST_OUTPUT_PATH#$PROJECT_DIR/}"
    OUTPUT_ENV="/workspace/project/${rel}"
  fi

  CONTAINER_PROGRAM="/workspace/project/assets/aggregate.codon"
  if [[ ! -f "${PROJECT_DIR}/assets/aggregate.codon" ]]; then
    CONTAINER_PROGRAM="/workspace/project/aggregate.codon"
  fi

  CONTAINER_NAME="syqure-${RUN_ID}-pid${PID}"
  docker rm -f "$CONTAINER_NAME" >/dev/null 2>&1 || true

  run_cmd=(docker run --rm --name "$CONTAINER_NAME" --platform "$PLATFORM"
    -e "SEQURE_TRANSPORT=file"
    -e "SEQURE_FILE_DIR=shared/syqure/${RUN_ID}"
    -e "SEQURE_FILE_POLL_MS=${POLL_MS}"
    -e "SEQURE_CP_COUNT=${CP_COUNT}"
    -e "SEQURE_PARTY_EMAILS=${PARTY_EMAILS}"
    -e "SEQURE_DATASITES_ROOT=${DATASITES_ROOT_IN_CONTAINER}"
    -e "SEQURE_LOCAL_EMAIL=${CURRENT_EMAIL}"
    ${SEQURE_FILE_KEEP:+-e "SEQURE_FILE_KEEP=${SEQURE_FILE_KEEP}"}
    ${SEQURE_FILE_DEBUG:+-e "SEQURE_FILE_DEBUG=${SEQURE_FILE_DEBUG}"}
    ${INPUT_ENV:+-e "BV_INPUT_MY_VAL=${INPUT_ENV}"}
    ${OUTPUT_ENV:+-e "BV_OUTPUT_AGGREGATE_RESULT=${OUTPUT_ENV}"}
    -v "${PROJECT_DIR}:/workspace/project"
    ${MOUNT_SPEC}
    "$IMAGE_NAME"
    syqure "$CONTAINER_PROGRAM" -- "$PID")

  if [[ -n "$LOG_PATH" ]]; then
    "${run_cmd[@]}" 2>&1 | tee "$LOG_PATH"
    exit "${PIPESTATUS[0]}"
  else
    exec "${run_cmd[@]}"
  fi
fi
