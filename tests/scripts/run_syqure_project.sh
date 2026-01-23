#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF' >&2
Usage: run_syqure_project.sh --sandbox DIR --bv PATH --agg MSG --client1 MSG --client2 MSG

Options:
  --sandbox DIR   Sandbox root containing datasite folders
  --bv PATH       Path to bv binary
  --agg MSG       Aggregator message ID
  --client1 MSG   Client1 message ID
  --client2 MSG   Client2 message ID
EOF
}

SANDBOX=""
BV_BIN=""
AGG_MSG=""
C1_MSG=""
C2_MSG=""
EXPECTED_CONTAINERS=3

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

cleanup_containers() {
  docker rm -f $(docker ps -aq --filter "name=syqure-") >/dev/null 2>&1 || true
}

start_tailer() {
  local email="$1"
  local label="$2"
  local start_ts="$3"
  local base="${SANDBOX}/${email}/datasites/${email}/shared/syqure"
  python3 - "$base" "$label" "$start_ts" <<'PY' &
import os
import sys
import time

base = sys.argv[1]
label = sys.argv[2]
start_ts = int(sys.argv[3])

deadline = time.time() + 50
log_path = None
while time.time() < deadline and log_path is None:
    if os.path.isdir(base):
        best = None
        best_mtime = 0
        for entry in os.listdir(base):
            candidate = os.path.join(base, entry, "file_transport.log")
            if not os.path.isfile(candidate):
                continue
            try:
                mtime = int(os.path.getmtime(candidate))
            except OSError:
                continue
            if mtime >= start_ts and mtime >= best_mtime:
                best = candidate
                best_mtime = mtime
        log_path = best
    if log_path is None:
        time.sleep(0.5)

if log_path is None:
    sys.exit(0)

with open(log_path, "r") as f:
    f.seek(0, os.SEEK_END)
    while True:
        line = f.readline()
        if not line:
            time.sleep(0.1)
            continue
        sys.stdout.write(f"[{label}] {line}")
        sys.stdout.flush()
PY
  echo $!
}

results_ready_party() {
  local email="$1"
  local msg="$2"
  local pid="$3"
  if [[ -f "${SANDBOX}/${email}/.biovault/runs/${msg}/results-test/${pid}_out.txt" ]]; then
    return 0
  fi
  compgen -G "${SANDBOX}/${email}/private/app_data/biovault/submissions/*/results-test/${pid}_out.txt" >/dev/null
}

results_ready() {
  results_ready_party "aggregator@sandbox.local" "$AGG_MSG" "0" &&
    results_ready_party "client1@sandbox.local" "$C1_MSG" "1" &&
    results_ready_party "client2@sandbox.local" "$C2_MSG" "2"
}

dump_container_logs() {
  echo "=== Syqure containers (all) ===" >&2
  docker ps -a --filter "name=syqure-" >&2 || true
  local ids
  ids=$(docker ps -aq --filter "name=syqure-") || true
  if [[ -n "$ids" ]]; then
    for id in $ids; do
      echo "--- docker logs $id ---" >&2
      docker logs "$id" >&2 || true
    done
  fi
}

dump_syqure_logs() {
  local email="$1"
  local msg="$2"
  local label="$3"
  local log_dir="${SANDBOX}/${email}/.biovault/runs/${msg}/syqure-logs"
  if [[ -d "$log_dir" ]]; then
    for log in "$log_dir"/*.log; do
      [[ -f "$log" ]] || continue
      echo "=== ${label} syqure log: $(basename "$log") ==="
      tail -n 200 "$log"
    done
  else
    echo "=== ${label} syqure logs missing: ${log_dir} ==="
  fi
}

cleanup_containers

start_ts="$(date +%s)"
tail_pids=()
tail_pids+=("$(start_tailer "aggregator@sandbox.local" "agg" "$start_ts")")
tail_pids+=("$(start_tailer "client1@sandbox.local" "c1" "$start_ts")")
tail_pids+=("$(start_tailer "client2@sandbox.local" "c2" "$start_ts")")

run_party() {
  local email="$1"
  local msg="$2"
  local home="${SANDBOX}/${email}"
  if [[ ! -d "$home" ]]; then
    echo "Missing datasite dir: ${home}" >&2
    return 1
  fi
  HOME="$home" \
    SYFTBOX_EMAIL="$email" \
    SYFTBOX_DATA_DIR="$home" \
    SYFTBOX_CONFIG_PATH="$home/.syftbox/config.json" \
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

for pid in "${tail_pids[@]}"; do
  kill "$pid" >/dev/null 2>&1 || true
done

if [[ "$status" -eq 0 ]] && ! results_ready; then
  echo "Syqure results not found after runs completed." >&2
  dump_container_logs
  exit 1
fi

dump_syqure_logs "aggregator@sandbox.local" "$AGG_MSG" "agg"
dump_syqure_logs "client1@sandbox.local" "$C1_MSG" "c1"
dump_syqure_logs "client2@sandbox.local" "$C2_MSG" "c2"

exit "$status"
