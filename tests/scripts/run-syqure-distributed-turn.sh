#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

TURN_CONTAINER="${TURN_CONTAINER:-biovault-turn-test}"
TURN_IMAGE="${TURN_IMAGE:-coturn/coturn:4.6.3}"
TURN_REALM="${TURN_REALM:-biovault.local}"
TURN_USER="${TURN_USER:-syftbox}"
TURN_PASS="${TURN_PASS:-syftbox}"
TURN_PORT="${TURN_PORT:-3478}"
TURN_HOST="${TURN_HOST:-127.0.0.1}"
KEEP_TURN="${KEEP_TURN:-0}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [--keep-turn]

Runs syqure distributed scenario in Docker with hotlink p2p-only and TURN configured.
This is a practical NAT/TURN-style validation (WebRTC-only, no WS fallback).

Options:
  --keep-turn   Do not remove TURN container after run
  -h, --help    Show this help

Env overrides:
  TURN_CONTAINER, TURN_IMAGE, TURN_REALM, TURN_USER, TURN_PASS, TURN_PORT, TURN_HOST, KEEP_TURN
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --keep-turn)
      KEEP_TURN=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
  shift
done

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

require_cmd docker

cleanup() {
  if [[ "${KEEP_TURN}" == "1" ]]; then
    echo "[turn] keeping container ${TURN_CONTAINER}"
    return
  fi
  docker rm -f "${TURN_CONTAINER}" >/dev/null 2>&1 || true
}
trap cleanup EXIT

if docker ps -a --format '{{.Names}}' | rg -x "${TURN_CONTAINER}" >/dev/null 2>&1; then
  echo "[turn] removing existing container ${TURN_CONTAINER}"
  docker rm -f "${TURN_CONTAINER}" >/dev/null
fi

echo "[turn] starting ${TURN_CONTAINER} (${TURN_IMAGE}) on ${TURN_HOST}:${TURN_PORT}"
docker run -d \
  --name "${TURN_CONTAINER}" \
  -p "${TURN_PORT}:${TURN_PORT}/udp" \
  -p "${TURN_PORT}:${TURN_PORT}/tcp" \
  "${TURN_IMAGE}" \
  -n \
  --realm="${TURN_REALM}" \
  --lt-cred-mech \
  --user="${TURN_USER}:${TURN_PASS}" \
  --listening-port="${TURN_PORT}" \
  --fingerprint \
  --no-cli \
  --log-file=stdout >/dev/null

echo "[turn] waiting for TURN logs"
for _ in $(seq 1 20); do
  if docker logs "${TURN_CONTAINER}" 2>&1 | rg -q "listening|startup complete|TURN Server"; then
    break
  fi
  sleep 1
done

cd "${ROOT_DIR}"

export SYFTBOX_HOTLINK_ICE_SERVERS="turn:${TURN_HOST}:${TURN_PORT}?transport=udp,turn:${TURN_HOST}:${TURN_PORT}?transport=tcp"
export SYFTBOX_HOTLINK_TURN_USER="${TURN_USER}"
export SYFTBOX_HOTLINK_TURN_PASS="${TURN_PASS}"
export BV_SYFTBOX_HOTLINK_ICE_SERVERS="${SYFTBOX_HOTLINK_ICE_SERVERS}"
export BV_SYFTBOX_HOTLINK_TURN_USER="${TURN_USER}"
export BV_SYFTBOX_HOTLINK_TURN_PASS="${TURN_PASS}"

echo "[run] ICE servers: ${SYFTBOX_HOTLINK_ICE_SERVERS}"
echo "[run] running syqure distributed (docker + hotlink p2p-only, ws fallback disabled)"

/usr/bin/time -p ./test-scenario.sh \
  --docker \
  --syqure-transport hotlink \
  --hotlink-p2p-only \
  tests/scenarios/syqure-distributed.yaml
