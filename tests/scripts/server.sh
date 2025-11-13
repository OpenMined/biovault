#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./server.sh [--reset]

Options:
  --reset    Stop SyftBox + MinIO containers and remove their volumes before starting.
EOF
}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DOCKER_DIR="$ROOT_DIR/syftbox/docker"
SERVER_URL="${SYFTBOX_SERVER_URL:-http://localhost:8080}"
RESET_STACK=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reset)
      RESET_STACK=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
  shift
done

command -v docker >/dev/null 2>&1 || { echo "docker is required" >&2; exit 1; }
command -v curl >/dev/null 2>&1 || { echo "curl is required" >&2; exit 1; }

[[ -d "$DOCKER_DIR" ]] || { echo "Missing SyftBox docker directory: $DOCKER_DIR" >&2; exit 1; }

if (( RESET_STACK )); then
  echo "Resetting SyftBox stack (containers + volumes)..."
  (cd "$DOCKER_DIR" && docker compose down -v || true)
fi

echo "Starting SyftBox server + MinIO..."
(cd "$DOCKER_DIR" && COMPOSE_BAKE=true docker compose up -d --build minio server)

echo "Waiting for $SERVER_URL to respond..."
for attempt in $(seq 1 60); do
  if curl -fsS --max-time 2 "$SERVER_URL" >/dev/null 2>&1; then
    echo "SyftBox server is ready."
    exit 0
  fi
  sleep 2
done

echo "Server did not respond within 120 seconds." >&2
exit 1
