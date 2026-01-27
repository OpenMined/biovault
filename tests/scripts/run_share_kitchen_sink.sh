#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF_USAGE' >&2
Usage: run_share_kitchen_sink.sh --sandbox DIR --bv PATH --agg MSG --client1 MSG --client2 MSG

Options:
  --sandbox DIR   Sandbox root containing datasite folders
  --bv PATH       Path to bv binary
  --agg MSG       Aggregator message ID
  --client1 MSG   Client1 message ID
  --client2 MSG   Client2 message ID
EOF_USAGE
}

SANDBOX=""
BV_BIN=""
AGG_MSG=""
C1_MSG=""
C2_MSG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sandbox) SANDBOX="${2:-}"; shift ;;
    --bv) BV_BIN="${2:-}"; shift ;;
    --agg) AGG_MSG="${2:-}"; shift ;;
    --client1) C1_MSG="${2:-}"; shift ;;
    --client2) C2_MSG="${2:-}"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
  shift
 done

if [[ -z "$SANDBOX" || -z "$BV_BIN" || -z "$AGG_MSG" || -z "$C1_MSG" || -z "$C2_MSG" ]]; then
  usage
  exit 1
fi

run_party() {
  local email="$1"
  local msg="$2"
  local home="${SANDBOX}/${email}"
  if [[ ! -d "$home" ]]; then
    echo "Missing datasite dir: ${home}" >&2
    return 1
  fi
  HOME="$home" \
    BIOVAULT_HOME="${home}/.biovault" \
    SYFTBOX_EMAIL="$email" \
    SYFTBOX_DATA_DIR="$home" \
    SYFTBOX_CONFIG_PATH="$home/.syftbox/config.json" \
    BV_BIN="$BV_BIN" \
    "$BV_BIN" message process "$msg" \
      --test \
      --approve \
      --non-interactive
}

pids=()
run_party "aggregator@sandbox.local" "$AGG_MSG" & pids+=("$!")
run_party "client1@sandbox.local" "$C1_MSG" & pids+=("$!")
run_party "client2@sandbox.local" "$C2_MSG" & pids+=("$!")

status=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    status=1
  fi
 done

results_ready_party() {
  local email="$1"
  local msg="$2"
  local pattern="$3"
  local base="${SANDBOX}/${email}/.biovault/runs/${msg}/results-test/pipeline"
  compgen -G "${base}/${pattern}" >/dev/null
}

results_ready() {
  results_ready_party "aggregator@sandbox.local" "$AGG_MSG" "collect/*/combined_hello.txt" && \
    results_ready_party "client1@sandbox.local" "$C1_MSG" "wait_clients/*/combined_from_agg.txt" && \
    results_ready_party "client2@sandbox.local" "$C2_MSG" "wait_clients/*/combined_from_agg.txt"
}

if [[ "$status" -eq 0 ]] && ! results_ready; then
  echo "Kitchen sink share results not found after runs completed." >&2
  exit 1
fi

exit "$status"
