#!/usr/bin/env bash
set -euo pipefail

# GitHub Actions Windows runners often provide `python` but not `python3` on PATH.
# Normalize so the rest of the script can keep using `python3`.
if ! command -v python3 >/dev/null 2>&1; then
  if command -v python >/dev/null 2>&1; then
    python3() { python "$@"; }
    export -f python3
  fi
fi

usage() {
  cat <<'EOF'
Usage: ./import-data.sh --csv FILE [--clients email1,email2] [--sandbox DIR] [--skip-detect] [--path-column NAME]

Options:
  --csv FILE         Samplesheet CSV to import (required).
  --clients list     Comma-separated list of sandbox client emails to target.
  --sandbox DIR      Override sandbox root (default: ./sandbox).
  --skip-detect      Skip 'bv files detect-csv' (assumes CSV already normalized).
  --path-column NAME Treat NAME as the path column in the source CSV (default: file_path).
  -h, --help         Show this message.

The script expects each client directory under the sandbox to be initialized via devstack.sh.
EOF
}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CLI_DIR="$ROOT_DIR/cli"
SANDBOX_DIR="${SANDBOX_DIR:-$ROOT_DIR/sandbox}"
CSV_PATH=""
RAW_CLIENTS=()
SKIP_DETECT=0
PATH_COLUMN="file_path"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --csv)
      CSV_PATH="${2:-}"
      [[ -z "$CSV_PATH" ]] && { echo "Missing value for --csv" >&2; usage >&2; exit 1; }
      shift
      ;;
    --clients)
      [[ $# -lt 2 ]] && { echo "Missing value for --clients" >&2; usage >&2; exit 1; }
      RAW_CLIENTS+=("$2")
      shift
      ;;
    --sandbox)
      [[ $# -lt 2 ]] && { echo "Missing value for --sandbox" >&2; usage >&2; exit 1; }
      SANDBOX_DIR="$2"
      shift
      ;;
    --skip-detect)
      SKIP_DETECT=1
      ;;
    --path-column)
      [[ $# -lt 2 ]] && { echo "Missing value for --path-column" >&2; usage >&2; exit 1; }
      PATH_COLUMN="$2"
      shift
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

require_bin() {
  # Use 'type' instead of 'command -v' to also find shell functions (e.g., python3 wrapper)
  type "$1" >/dev/null 2>&1 || { echo "Missing required tool: $1" >&2; exit 1; }
}

abs_path() {
  python3 - <<'PY' "$1"
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
}

[[ -n "$CSV_PATH" ]] || { echo "--csv is required" >&2; usage >&2; exit 1; }
require_bin python3
CSV_PATH="$(abs_path "$CSV_PATH")"
[[ -f "$CSV_PATH" ]] || { echo "CSV not found: $CSV_PATH" >&2; exit 1; }
SOURCE_CSV="$CSV_PATH"

require_bin cargo
[[ -d "$CLI_DIR" ]] || { echo "Missing cli checkout at $CLI_DIR" >&2; exit 1; }
[[ -d "$SANDBOX_DIR" ]] || { echo "Sandbox directory not found: $SANDBOX_DIR" >&2; exit 1; }

declare -a CLIENTS=()
add_client() {
  local raw trimmed lowered
  raw="$1"
  trimmed="${raw#"${raw%%[![:space:]]*}"}"
  trimmed="${trimmed%"${trimmed##*[![:space:]]}"}"
  [[ -z "$trimmed" ]] && return
  lowered="$(printf '%s' "$trimmed" | tr '[:upper:]' '[:lower:]')"
  if ((${#CLIENTS[@]})); then
    for existing in "${CLIENTS[@]}"; do
      if [[ "$existing" == "$lowered" ]]; then
        return
      fi
    done
  fi
  CLIENTS+=("$lowered")
}

if ((${#RAW_CLIENTS[@]})); then
  for block in "${RAW_CLIENTS[@]}"; do
    IFS=',' read -r -a pieces <<< "$block"
    for piece in "${pieces[@]}"; do
      add_client "$piece"
    done
  done
fi

if ((${#CLIENTS[@]} == 0)); then
  while IFS= read -r dir; do
    add_client "$(basename "$dir")"
  done < <(find "$SANDBOX_DIR" -mindepth 1 -maxdepth 1 -type d -print 2>/dev/null | sort)
fi

if ((${#CLIENTS[@]} == 0)); then
  echo "No sandbox clients found under $SANDBOX_DIR" >&2
  exit 1
fi

ensure_bv_binary() {
  local target="$CLI_DIR/target/release/bv"
  if [[ ! -x "$target" ]]; then
    echo "Building BioVault CLI..."
    (cd "$CLI_DIR" && cargo build --release >/dev/null)
  fi
  printf '%s\n' "$target"
}

BV_BIN="$(ensure_bv_binary)"

run_bv() {
  local client_dir="$1"; shift
  local email
  email="$(basename "$client_dir")"
  local data_dir="$client_dir"
  local config_path="$client_dir/.syftbox/config.json"
  local config_yaml="$client_dir/config.yaml"
  [[ -f "$config_path" ]] || { echo "Missing SyftBox config for $email at $config_path" >&2; exit 1; }
  HOME="$client_dir" \
  BIOVAULT_HOME="$([[ -f "$config_yaml" ]] && printf '%s' "$client_dir" || printf '%s' "$client_dir/.biovault")" \
  SYFTBOX_EMAIL="$email" \
  SYFTBOX_DATA_DIR="$data_dir" \
  SYFTBOX_CONFIG_PATH="$config_path" \
  "$BV_BIN" "$@"
}

count_csv_rows() {
  python3 - <<'PY' "$1"
import csv, sys
path = sys.argv[1]
with open(path, newline='') as fh:
    reader = csv.reader(fh)
    next(reader, None)
    rows = sum(1 for row in reader if any(col.strip() for col in row))
print(rows)
PY
}

prepare_csv_for_client() {
  local src="$1" dest="$2"
  python3 - <<'PY' "$src" "$dest"
import csv, sys, os
src, dest = sys.argv[1], sys.argv[2]
target = "file_path"
aliases = ["file_path", "genotype_file", "genotype_file_path", "path"]
src_dir = os.path.dirname(os.path.abspath(src))

with open(src, newline='', encoding='utf-8') as infile:
    reader = csv.reader(infile)
    rows = list(reader)
if not rows:
    print(f"{src} is empty", file=sys.stderr)
    sys.exit(1)
header = rows[0]
if target not in header:
    replacement = None
    for alias in aliases:
        if alias in header:
            replacement = alias
            break
    if replacement is None:
        print(f"CSV header missing any of {aliases}", file=sys.stderr)
        sys.exit(1)
    header = [target if col == replacement else col for col in header]
    path_col_idx = header.index(target)
else:
    path_col_idx = header.index(target)

rows[0] = header

# Resolve relative paths to absolute paths based on CSV directory
for i in range(1, len(rows)):
    if path_col_idx < len(rows[i]):
        path = rows[i][path_col_idx]
        if path and not os.path.isabs(path):
            # Resolve relative to source CSV directory
            rows[i][path_col_idx] = os.path.abspath(os.path.join(src_dir, path))

with open(dest, 'w', newline='', encoding='utf-8') as outfile:
    writer = csv.writer(outfile)
    writer.writerows(rows)
PY
}

extract_total_from_json() {
  python3 - <<'PY'
import json, sys
payload = sys.stdin.read()
start = None
for idx, ch in enumerate(payload):
    if ch in '{[':
        start = idx
        break
if start is None:
    print(0)
    sys.exit(1)
data = json.loads(payload[start:])
if isinstance(data, dict) and 'data' in data and isinstance(data['data'], dict):
    data = data['data']
files = data.get('files')
if isinstance(files, list):
    total = data.get('total', len(files))
else:
    total = data.get('total', 0)
print(int(total))
PY
}

EXPECTED_ROWS="$(count_csv_rows "$SOURCE_CSV")"
if [[ "$EXPECTED_ROWS" -lt 1 ]]; then
  echo "CSV $CSV_PATH contains no data rows after header" >&2
  exit 1
fi

echo "Importing data from $CSV_PATH into ${#CLIENTS[@]} client(s) (expected rows: $EXPECTED_ROWS)"

for email in "${CLIENTS[@]}"; do
  client_dir="$SANDBOX_DIR/$email"
  if [[ ! -d "$client_dir" ]]; then
    echo "Skipping $email (directory not found at $client_dir)" >&2
    continue
  fi

  echo ""
  echo "=== Client: $email ==="
  client_tmp_dir="$client_dir/.sandbox-import"
  mkdir -p "$client_tmp_dir"
  client_csv="$client_tmp_dir/$(basename "$SOURCE_CSV")"
  if ! prepare_csv_for_client "$SOURCE_CSV" "$client_csv"; then
    echo "Failed to prepare CSV for $email" >&2
    exit 1
  fi

  if (( SKIP_DETECT == 0 )); then
    if detect_out=$(run_bv "$client_dir" files detect-csv "$client_csv" -o "$client_csv" 2>&1); then
      echo "detect-csv:"
      echo "$detect_out"
    else
      echo "detect-csv failed for $email:"
      echo "$detect_out"
      exit 1
    fi
  else
    echo "Skipping detect-csv for $email"
  fi

  if import_out=$(run_bv "$client_dir" files import-csv "$client_csv" --non-interactive --format json 2>&1); then
    echo "import-csv:"
    echo "$import_out"
  else
    echo "import-csv failed for $email:"
    echo "$import_out"
    exit 1
  fi

  if list_out=$(run_bv "$client_dir" files list --format json 2>&1); then
    total="$(printf '%s' "$list_out" | extract_total_from_json 2>/dev/null || echo 0)"
    [[ "$total" =~ ^[0-9]+$ ]] || total=0
    echo "files list total: $total"
    if (( total < EXPECTED_ROWS )); then
      echo "⚠ Warning: catalog total ($total) is less than expected rows ($EXPECTED_ROWS)" >&2
    else
      echo "✓ Catalog contains at least the expected number of rows."
    fi
  else
    echo "files list failed for $email:"
    echo "$list_out"
    exit 1
  fi

  rm -f "$client_csv"

done

echo ""
echo "Import complete."
