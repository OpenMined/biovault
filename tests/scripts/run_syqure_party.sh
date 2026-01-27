#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF' >&2
Usage: run_syqure_party.sh --pid N --email EMAIL --run-id ID --party-emails CSV --datasites-root DIR [options]

Options:
  --image NAME        Docker image (default: syqure-cli-files)
  --platform PLAT     Docker platform (default: linux/amd64)
  --program PATH      Program path in container (default: /workspace/example/two_party_sum_tcp.codon)
  --example-dir DIR   Host example dir to mount (default: <repo>/syqure/example)
  --name NAME         Container name (default: syqure-<email>)
  --poll-ms MS        File poll interval (default: 50)
  --transport MODE    Transport (default: file)
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Detect Windows (Git Bash / MSYS / Cygwin)
IS_WINDOWS=0
if [[ "$(uname -s)" == MINGW* ]] || [[ "$(uname -s)" == MSYS* ]] || [[ "$(uname -s)" == CYGWIN* ]]; then
  IS_WINDOWS=1
fi

docker_host_path() {
  local path="$1"
  if (( IS_WINDOWS )); then
    if command -v cygpath >/dev/null 2>&1; then
      cygpath -m "$path"
      return
    fi
  fi
  echo "$path"
}

PID=""
EMAIL=""
RUN_ID=""
PARTY_EMAILS=""
DATASITES_ROOT=""
IMAGE_NAME="${IMAGE_NAME:-syqure-cli-files}"
PLATFORM="${PLATFORM:-linux/amd64}"
PROGRAM_PATH="${PROGRAM_PATH:-/workspace/example/two_party_sum_tcp.codon}"
EXAMPLE_DIR="${EXAMPLE_DIR:-$ROOT_DIR/../syqure/example}"
CONTAINER_NAME=""
POLL_MS="${POLL_MS:-50}"
TRANSPORT="${TRANSPORT:-file}"
FILE_KEEP="${SEQURE_FILE_KEEP:-}"
FILE_DEBUG="${SEQURE_FILE_DEBUG:-}"
CP_COUNT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pid) PID="${2:-}"; shift ;;
    --email) EMAIL="${2:-}"; shift ;;
    --run-id) RUN_ID="${2:-}"; shift ;;
    --party-emails) PARTY_EMAILS="${2:-}"; shift ;;
    --datasites-root) DATASITES_ROOT="${2:-}"; shift ;;
    --image) IMAGE_NAME="${2:-}"; shift ;;
    --platform) PLATFORM="${2:-}"; shift ;;
    --program) PROGRAM_PATH="${2:-}"; shift ;;
    --example-dir) EXAMPLE_DIR="${2:-}"; shift ;;
    --name) CONTAINER_NAME="${2:-}"; shift ;;
    --poll-ms) POLL_MS="${2:-}"; shift ;;
    --transport) TRANSPORT="${2:-}"; shift ;;
    --cp-count) CP_COUNT="${2:-}"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
  shift
done

if [[ -z "$PID" || -z "$EMAIL" || -z "$RUN_ID" || -z "$PARTY_EMAILS" || -z "$DATASITES_ROOT" ]]; then
  usage
  exit 1
fi

if [[ -z "$CP_COUNT" ]]; then
  IFS=',' read -r -a parties <<< "$PARTY_EMAILS"
  CP_COUNT="${#parties[@]}"
fi

if [[ -z "$CONTAINER_NAME" ]]; then
  safe_email="${EMAIL//@/-}"
  safe_email="${safe_email//./-}"
  CONTAINER_NAME="syqure-${safe_email}-${RUN_ID}"
fi

mkdir -p "$DATASITES_ROOT"
docker rm -f "$CONTAINER_NAME" >/dev/null 2>&1 || true

waited=0
while [[ ! -d "$DATASITES_ROOT" && $waited -lt 50 ]]; do
  sleep 0.1
  waited=$((waited + 1))
done

DATASITES_ROOT_IN_CONTAINER="/datasites"
DATASITES_ROOT_HOST="$(docker_host_path "$DATASITES_ROOT")"
MOUNT_SPEC="-v ${DATASITES_ROOT_HOST}:${DATASITES_ROOT_IN_CONTAINER}"
if ! docker run --rm ${MOUNT_SPEC} alpine:3.19 sh -c 'test -d /datasites' >/dev/null 2>&1; then
  parent_dir="$(dirname "$DATASITES_ROOT")"
  if [[ ! -d "$parent_dir" ]]; then
    echo "Datasites parent directory not found: $parent_dir" >&2
    exit 1
  fi
  parent_dir_host="$(docker_host_path "$parent_dir")"
  DATASITES_ROOT_IN_CONTAINER="/datasites-root/datasites"
  MOUNT_SPEC="-v ${parent_dir_host}:/datasites-root"
fi

if [[ "$PROGRAM_PATH" == /workspace/example/* ]]; then
  rel="${PROGRAM_PATH#/workspace/example/}"
  if [ ! -f "${EXAMPLE_DIR}/${rel}" ]; then
    echo "Example file not found at ${EXAMPLE_DIR}/${rel}" >&2
    exit 1
  fi
fi

if [[ -n "$FILE_DEBUG" ]]; then
  log_path="${DATASITES_ROOT_IN_CONTAINER}/${EMAIL}/shared/syqure/${RUN_ID}/file_transport.log"
  example_dir_host="$(docker_host_path "$EXAMPLE_DIR")"
  if (( IS_WINDOWS )); then
    MSYS2_ARG_CONV_EXCL="*" docker run -d --name "$CONTAINER_NAME" --platform "$PLATFORM" \
      -e "SEQURE_TRANSPORT=${TRANSPORT}" \
      -e "SEQURE_FILE_DIR=shared/syqure/${RUN_ID}" \
      -e "SEQURE_FILE_POLL_MS=${POLL_MS}" \
      -e "SEQURE_CP_COUNT=${CP_COUNT}" \
      -e "SEQURE_PARTY_EMAILS=${PARTY_EMAILS}" \
      -e "SEQURE_DATASITES_ROOT=${DATASITES_ROOT_IN_CONTAINER}" \
      -e "SEQURE_LOCAL_EMAIL=${EMAIL}" \
      ${FILE_KEEP:+-e "SEQURE_FILE_KEEP=${FILE_KEEP}"} \
      ${FILE_DEBUG:+-e "SEQURE_FILE_DEBUG=${FILE_DEBUG}"} \
      -v "${example_dir_host}:/workspace/example:ro" \
      ${MOUNT_SPEC} \
      "$IMAGE_NAME" sh -c "syqure \"$PROGRAM_PATH\" -- \"$PID\" & pid=\$!; \
        log=\"${log_path}\"; \
        i=0; while [ \$i -lt 300 ]; do [ -f \"\$log\" ] && break; sleep 0.2; i=\$((i+1)); done; \
        if [ -f \"\$log\" ]; then tail -n 50 -F \"\$log\" & tpid=\$!; fi; \
        wait \$pid; status=\$?; \
        if [ -n \"\${tpid:-}\" ]; then kill \$tpid 2>/dev/null || true; fi; \
        exit \$status"
  else
    docker run -d --name "$CONTAINER_NAME" --platform "$PLATFORM" \
      -e "SEQURE_TRANSPORT=${TRANSPORT}" \
      -e "SEQURE_FILE_DIR=shared/syqure/${RUN_ID}" \
      -e "SEQURE_FILE_POLL_MS=${POLL_MS}" \
      -e "SEQURE_CP_COUNT=${CP_COUNT}" \
      -e "SEQURE_PARTY_EMAILS=${PARTY_EMAILS}" \
      -e "SEQURE_DATASITES_ROOT=${DATASITES_ROOT_IN_CONTAINER}" \
      -e "SEQURE_LOCAL_EMAIL=${EMAIL}" \
      ${FILE_KEEP:+-e "SEQURE_FILE_KEEP=${FILE_KEEP}"} \
      ${FILE_DEBUG:+-e "SEQURE_FILE_DEBUG=${FILE_DEBUG}"} \
      -v "${EXAMPLE_DIR}:/workspace/example:ro" \
      ${MOUNT_SPEC} \
      "$IMAGE_NAME" sh -c "syqure \"$PROGRAM_PATH\" -- \"$PID\" & pid=\$!; \
        log=\"${log_path}\"; \
        i=0; while [ \$i -lt 300 ]; do [ -f \"\$log\" ] && break; sleep 0.2; i=\$((i+1)); done; \
        if [ -f \"\$log\" ]; then tail -n 50 -F \"\$log\" & tpid=\$!; fi; \
        wait \$pid; status=\$?; \
        if [ -n \"\${tpid:-}\" ]; then kill \$tpid 2>/dev/null || true; fi; \
        exit \$status"
  fi
  exit $?
fi

example_dir_host="$(docker_host_path "$EXAMPLE_DIR")"
if (( IS_WINDOWS )); then
  MSYS2_ARG_CONV_EXCL="*" docker run -d --name "$CONTAINER_NAME" --platform "$PLATFORM" \
    -e "SEQURE_TRANSPORT=${TRANSPORT}" \
    -e "SEQURE_FILE_DIR=shared/syqure/${RUN_ID}" \
    -e "SEQURE_FILE_POLL_MS=${POLL_MS}" \
    -e "SEQURE_CP_COUNT=${CP_COUNT}" \
    -e "SEQURE_PARTY_EMAILS=${PARTY_EMAILS}" \
    -e "SEQURE_DATASITES_ROOT=${DATASITES_ROOT_IN_CONTAINER}" \
    ${FILE_KEEP:+-e "SEQURE_FILE_KEEP=${FILE_KEEP}"} \
    ${FILE_DEBUG:+-e "SEQURE_FILE_DEBUG=${FILE_DEBUG}"} \
    -v "${example_dir_host}:/workspace/example:ro" \
    ${MOUNT_SPEC} \
    "$IMAGE_NAME" syqure "$PROGRAM_PATH" -- "$PID"
else
  docker run -d --name "$CONTAINER_NAME" --platform "$PLATFORM" \
    -e "SEQURE_TRANSPORT=${TRANSPORT}" \
    -e "SEQURE_FILE_DIR=shared/syqure/${RUN_ID}" \
    -e "SEQURE_FILE_POLL_MS=${POLL_MS}" \
    -e "SEQURE_CP_COUNT=${CP_COUNT}" \
    -e "SEQURE_PARTY_EMAILS=${PARTY_EMAILS}" \
    -e "SEQURE_DATASITES_ROOT=${DATASITES_ROOT_IN_CONTAINER}" \
    ${FILE_KEEP:+-e "SEQURE_FILE_KEEP=${FILE_KEEP}"} \
    ${FILE_DEBUG:+-e "SEQURE_FILE_DEBUG=${FILE_DEBUG}"} \
    -v "${EXAMPLE_DIR}:/workspace/example:ro" \
    ${MOUNT_SPEC} \
    "$IMAGE_NAME" syqure "$PROGRAM_PATH" -- "$PID"
fi
