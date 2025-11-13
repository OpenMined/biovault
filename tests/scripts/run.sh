#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./run.sh --datasite EMAIL --project PATH [options]

Options:
  --datasite EMAIL     Sandbox client email to run as (required).
  --project PATH       Pipeline or project directory/file (required).
  --sandbox DIR        Override sandbox root (default: ./sandbox).
  --samplesheet FILE   Use an existing samplesheet instead of generating one.
  --data VALUE         Select participants when generating samplesheet ('all' or comma list).
  --results DIR        Directory for pipeline results (default: .sandbox-run/results/<timestamp>).
  --set KEY=VALUE      Pass through additional --set inputs.* overrides (repeatable).
  --test               Run pipeline with --test flag.
  --dry-run            Pass --dry-run to bv run.
  --with-docker        Pass --with-docker to bv run.
  --skip-generate      Skip automatic samplesheet generation (requires --samplesheet).
  -h, --help           Show this message.

The script activates the specified sandbox client via sbenv, generates a samplesheet
covering all imported files for that datasite (unless --samplesheet is provided),
and runs `bv run` for the supplied project/pipeline with the samplesheet wired to
`inputs.samplesheet`.
EOF
}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CLI_DIR="$ROOT_DIR/cli"
SBENV_BIN="$ROOT_DIR/sbenv/sbenv"
SANDBOX_DIR="${SANDBOX_DIR:-$ROOT_DIR/sandbox}"

DATASITE=""
PROJECT_ARG=""
SAMPLESHEET_OVERRIDE=""
RESULTS_DIR=""
SKIP_GENERATE=0
RUN_TEST=0
RUN_DRY=0
RUN_DOCKER=0
declare -a EXTRA_SET_ARGS=()
DATA_FILTER="all"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --datasite)
      DATASITE="${2:-}"
      [[ -z "$DATASITE" ]] && { echo "Missing value for --datasite" >&2; usage >&2; exit 1; }
      shift
      ;;
    --project)
      PROJECT_ARG="${2:-}"
      [[ -z "$PROJECT_ARG" ]] && { echo "Missing value for --project" >&2; usage >&2; exit 1; }
      shift
      ;;
    --sandbox)
      SANDBOX_DIR="${2:-}"
      [[ -z "$SANDBOX_DIR" ]] && { echo "Missing value for --sandbox" >&2; usage >&2; exit 1; }
      shift
      ;;
    --samplesheet)
      SAMPLESHEET_OVERRIDE="${2:-}"
      [[ -z "$SAMPLESHEET_OVERRIDE" ]] && { echo "Missing value for --samplesheet" >&2; usage >&2; exit 1; }
      shift
      ;;
    --data)
      DATA_FILTER="${2:-}"
      [[ -z "$DATA_FILTER" ]] && { echo "Missing value for --data" >&2; usage >&2; exit 1; }
      shift
      ;;
    --results)
      RESULTS_DIR="${2:-}"
      [[ -z "$RESULTS_DIR" ]] && { echo "Missing value for --results" >&2; usage >&2; exit 1; }
      shift
      ;;
    --set)
      [[ $# -lt 2 ]] && { echo "Missing value for --set" >&2; usage >&2; exit 1; }
      EXTRA_SET_ARGS+=("$2")
      shift
      ;;
    --test)
      RUN_TEST=1
      ;;
    --dry-run)
      RUN_DRY=1
      ;;
    --with-docker)
      RUN_DOCKER=1
      ;;
    --skip-generate)
      SKIP_GENERATE=1
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
  command -v "$1" >/dev/null 2>&1 || { echo "Missing required tool: $1" >&2; exit 1; }
}

abs_path() {
  python3 - <<'PY' "$1"
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
}

[[ -n "$DATASITE" ]] || { echo "--datasite is required" >&2; usage >&2; exit 1; }
[[ -n "$PROJECT_ARG" ]] || { echo "--project is required" >&2; usage >&2; exit 1; }

require_bin python3
require_bin cargo
[[ -x "$SBENV_BIN" ]] || { echo "Missing sbenv binary at $SBENV_BIN" >&2; exit 1; }
[[ -d "$CLI_DIR" ]] || { echo "Missing cli directory at $CLI_DIR" >&2; exit 1; }

SANDBOX_DIR="$(abs_path "$SANDBOX_DIR")"
CLIENT_DIR="$SANDBOX_DIR/$DATASITE"
[[ -d "$CLIENT_DIR" ]] || { echo "Sandbox client not found at $CLIENT_DIR" >&2; exit 1; }

PROJECT_PATH="$(abs_path "$PROJECT_ARG")"
if [[ -d "$PROJECT_PATH" ]]; then
  :
elif [[ -f "$PROJECT_PATH" ]]; then
  :
else
  echo "Project path not found: $PROJECT_PATH" >&2
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
  (
    set +u
    cd "$client_dir"
    eval "$("$SBENV_BIN" activate --quiet)"
    "$BV_BIN" "$@"
  )
}

timestamp() {
  date -u +%Y%m%dT%H%M%SZ
}

if [[ -n "$SAMPLESHEET_OVERRIDE" ]]; then
  SAMPLESHEET_PATH="$(abs_path "$SAMPLESHEET_OVERRIDE")"
  [[ -f "$SAMPLESHEET_PATH" ]] || { echo "Samplesheet not found: $SAMPLESHEET_PATH" >&2; exit 1; }
else
  if (( SKIP_GENERATE )); then
    echo "--skip-generate requires --samplesheet" >&2
    exit 1
  fi
  SAMPLE_DIR="$CLIENT_DIR/.sandbox-run"
  mkdir -p "$SAMPLE_DIR"
  SAMPLESHEET_PATH="$SAMPLE_DIR/samplesheet-$(timestamp).csv"
fi

generate_samplesheet() {
  local output="$1"
  local filter_expr="$2"
  local lower
  lower="$(printf '%s' "$filter_expr" | tr '[:upper:]' '[:lower:]')"
  local args=("samplesheet" "export-catalog" "$output" "--path-column" "genotype_file" "--data-type" "genotype")
  if [[ "$lower" != "all" && -n "$filter_expr" ]]; then
    args+=("--participants" "$filter_expr")
  fi
  if ! run_bv "$CLIENT_DIR" "${args[@]}"; then
    echo "Failed to export samplesheet for $DATASITE" >&2
    exit 1
  fi
}

if [[ -z "$SAMPLESHEET_OVERRIDE" ]]; then
  generate_samplesheet "$SAMPLESHEET_PATH" "$DATA_FILTER"
fi

echo "[DEBUG] Creating results directory..."
if [[ -z "$RESULTS_DIR" ]]; then
  RESULTS_DIR="$CLIENT_DIR/.sandbox-run/results/$(basename "$PROJECT_PATH")-$(timestamp)"
fi
RESULTS_DIR="$(abs_path "$RESULTS_DIR")"
echo "[DEBUG] Results dir: $RESULTS_DIR"
mkdir -p "$RESULTS_DIR"
echo "[DEBUG] Results directory created"

infer_data_dir() {
  local sheet="$1"
  echo "[DEBUG infer_data_dir] Processing: $sheet" >&2
  python3 - <<'PY' "$sheet"
import csv, os, sys
sheet = sys.argv[1]
with open(sheet, newline='', encoding='utf-8') as fh:
    reader = csv.reader(fh)
    try:
        header = next(reader)
    except StopIteration:
        print("Provided samplesheet is empty", file=sys.stderr)
        sys.exit(1)
columns = [col.strip() for col in header]
path_col = None
for candidate in ("file_path", "genotype_file", "genotype_file_path", "path"):
    if candidate in columns:
        path_col = candidate
        break
if path_col is None:
    print("Samplesheet missing recognizable path column", file=sys.stderr)
    sys.exit(1)
idx = columns.index(path_col)
paths = []
with open(sheet, newline='', encoding='utf-8') as fh:
    reader = csv.reader(fh)
    next(reader, None)
    for row in reader:
        if idx >= len(row):
            continue
        value = row[idx].strip()
        if not value:
            continue
        abs_path = os.path.abspath(value)
        paths.append(abs_path)
if not paths:
    print("Samplesheet contains no valid file paths", file=sys.stderr)
    sys.exit(1)
common = os.path.commonpath(paths)
print(common)
PY
}

requires_data_dir() {
  local target="$1"
  echo "[DEBUG requires_data_dir] Checking: $target" >&2
  local file=""
  if [[ -f "$target" ]]; then
    file="$target"
    echo "[DEBUG requires_data_dir] Found file: $file" >&2
  elif [[ -d "$target" ]]; then
    echo "[DEBUG requires_data_dir] Target is directory" >&2
    if [[ -f "$target/pipeline.yaml" ]]; then
      file="$target/pipeline.yaml"
      echo "[DEBUG requires_data_dir] Found: pipeline.yaml" >&2
    elif [[ -f "$target/project.yaml" ]]; then
      file="$target/project.yaml"
      echo "[DEBUG requires_data_dir] Found: project.yaml" >&2
    fi
  fi
  if [[ -z "$file" ]]; then
    echo "[DEBUG requires_data_dir] No file found, returning false" >&2
    return 1
  fi
  echo "[DEBUG requires_data_dir] Searching for data_dir in: $file" >&2
  if command -v rg >/dev/null 2>&1; then
    rg -q --fixed-strings --word-regexp "data_dir" "$file"
  else
    grep -q "data_dir" "$file"
  fi
  local result=$?
  echo "[DEBUG requires_data_dir] Search result: $result" >&2
  return $result
}

echo "[DEBUG] Checking samplesheet arg..."
SAMPLE_SET_PRESENT=0
if ((${#EXTRA_SET_ARGS[@]:-0})); then
  for entry in "${EXTRA_SET_ARGS[@]}"; do
    if [[ "$entry" == inputs.samplesheet=* ]]; then
      SAMPLE_SET_PRESENT=1
      break
    fi
  done
fi
if (( SAMPLE_SET_PRESENT == 0 )); then
  EXTRA_SET_ARGS+=("inputs.samplesheet=$SAMPLESHEET_PATH")
fi
echo "[DEBUG] Samplesheet arg set: $SAMPLESHEET_PATH"

echo "[DEBUG] Checking data_dir arg..."
DATA_SET_PRESENT=0
if ((${#EXTRA_SET_ARGS[@]:-0})); then
  for entry in "${EXTRA_SET_ARGS[@]}"; do
    if [[ "$entry" == inputs.data_dir=* ]]; then
      DATA_SET_PRESENT=1
      break
    fi
  done
fi
echo "[DEBUG] Data dir present: $DATA_SET_PRESENT"
echo "[DEBUG] Checking if project requires data_dir..."
if (( DATA_SET_PRESENT == 0 )) && requires_data_dir "$PROJECT_PATH"; then
  echo "[DEBUG] Project requires data_dir, inferring from samplesheet..."
  if inferred_dir="$(infer_data_dir "$SAMPLESHEET_PATH")"; then
    EXTRA_SET_ARGS+=("inputs.data_dir=$inferred_dir")
    echo "Inferred data dir: $inferred_dir"
  else
    echo "Failed to infer data directory from $SAMPLESHEET_PATH" >&2
    exit 1
  fi
else
  echo "[DEBUG] Project does not require data_dir or already set"
fi

echo "[DEBUG] Building run command..."
RUN_CMD=(run "$PROJECT_PATH" --results-dir "$RESULTS_DIR")
(( RUN_TEST )) && RUN_CMD+=(--test)
(( RUN_DRY )) && RUN_CMD+=(--dry-run)
(( RUN_DOCKER )) && RUN_CMD+=(--with-docker)
for entry in "${EXTRA_SET_ARGS[@]}"; do
  RUN_CMD+=(--set "$entry")
done
echo "[DEBUG] Run command built successfully"

echo "Running pipeline:"
printf '  bv %s\n' "${RUN_CMD[*]}"

if ! run_bv "$CLIENT_DIR" "${RUN_CMD[@]}"; then
  echo "Pipeline run failed. Check logs under $RESULTS_DIR" >&2
  exit 1
fi

echo ""
echo "Pipeline completed successfully."
echo "  Datasite : $DATASITE"
echo "  Project  : $PROJECT_PATH"
echo "  Results  : $RESULTS_DIR"
echo "  Samples  : $SAMPLESHEET_PATH"
