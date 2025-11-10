#!/usr/bin/env bash
set -euo pipefail

WAIT_VALUE="false"

ensure_submodules_populated() {
  if [[ ! -f .gitmodules ]]; then
    return
  fi

  # Collect submodule paths defined in .gitmodules
  local submodules=()
  while IFS=' ' read -r _ path; do
    [[ -n "$path" ]] && submodules+=("$path")
  done < <(git config --file .gitmodules --get-regexp path 2>/dev/null || true)

  if [[ ${#submodules[@]} -eq 0 ]]; then
    return
  fi

  local missing=()
  for path in "${submodules[@]}"; do
    if [[ ! -d "$path" ]]; then
      missing+=("$path (directory missing)")
      continue
    fi
    if ! find "$path" -mindepth 1 -maxdepth 1 -not -name '.git' -print -quit | grep -q .; then
      missing+=("$path (empty)")
    fi
  done

  if ((${#missing[@]})); then
    echo "The following submodule directories appear empty:" >&2
    for entry in "${missing[@]}"; do
      echo "  - $entry" >&2
    done
    echo "Run ./submodules.sh to fetch their contents before running tests." >&2
    exit 1
  fi
}

while (($#)); do
  case "$1" in
    --wait)
      WAIT_VALUE="1"
      if (($# > 1)) && [[ "$2" =~ ^[0-9]+$ ]]; then
        WAIT_VALUE="$2"
        shift
      fi
      ;;
    --wait=*)
      value="${1#*=}"
      if [[ -z "$value" ]]; then
        WAIT_VALUE="1"
      elif [[ "$value" =~ ^[0-9]+$ ]]; then
        WAIT_VALUE="$value"
      else
        echo "Invalid value for --wait: $value" >&2
        exit 1
      fi
      ;;
    -h|--help)
      cat <<'USAGE'
Usage: ./test-integration.sh [--wait[=STEP]]

Options:
  --wait[=STEP]    Pause between integration steps starting at STEP (default STEP=1).
USAGE
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -* )
      echo "Unknown option: $1" >&2
      exit 1
      ;;
    *)
      echo "Unknown positional argument: $1" >&2
      exit 1
      ;;
  esac
  shift
done

ensure_submodules_populated

rm -rf test-clients-local

if command -v docker >/dev/null 2>&1; then
  # Stop and remove only the containers that the integration test uses
  containers=(
    "syftbox-server"
    "syftbox-minio"
  )

  for name in "${containers[@]}"; do
    if docker ps -a --format '{{.Names}}' | grep -Fxq "$name"; then
      docker rm --force "$name" || true
    fi
  done

  # Remove the named volume if present
  if docker volume inspect docker_minio-data >/dev/null 2>&1; then
    docker volume rm docker_minio-data || true
  fi
else
  echo "Docker not available; skipping container cleanup" >&2
fi

# Run the integration loop
if [[ "$WAIT_VALUE" != "false" ]]; then
  just test-integration-local-inspect "$WAIT_VALUE"
else
  just test-integration-local-inspect
fi
