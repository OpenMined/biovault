#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${ROOT_DIR}/licenses"

mkdir -p "${OUT_DIR}"
export OUT_DIR

IGNORE_DIRS=("syftbox" "syftbox-sdk" "biovault" "beaver" "sbenv" "bioscript")
IGNORE_DIRS_CSV="$(IFS=','; echo "${IGNORE_DIRS[*]}")"
export IGNORE_DIRS_CSV

INSTALL_MISSING=false
SCAN_TARGETS="all"
FETCH_LICENSES=false

usage() {
  cat <<'EOF'
Usage: ./licenses.sh [--install] [--scan all|go|rust|go,rust] [--fetch]

Options:
  --install   Install missing license tools for detected languages
  --scan      Restrict scans to a subset (default: all)
  --fetch     Fetch and cache license texts for Go modules (requires network)
EOF
}

find_manifests() {
  local pattern="$1"
  local -a prune_expr=()

  for dir in "${IGNORE_DIRS[@]}"; do
    prune_expr+=(-path "${ROOT_DIR}/${dir}" -o -path "${ROOT_DIR}/${dir}/*" -o)
  done

  if (( ${#prune_expr[@]} > 0 )); then
    prune_expr+=(-false)
    find "${ROOT_DIR}" \( "${prune_expr[@]}" \) -prune -o -name "${pattern}" -print
  else
    find "${ROOT_DIR}" -name "${pattern}" -print
  fi
}

slugify_path() {
  local path="$1"
  local rel="${path#${ROOT_DIR}/}"
  if [[ "${rel}" == "${path}" ]]; then
    rel="${path}"
  fi
  if [[ -z "${rel}" ]]; then
    rel="root"
  fi
  rel="${rel%/}"
  rel="${rel//\//_}"
  rel="${rel//./_}"
  echo "${rel}"
}

filter_rust_licenses() {
  local rust_file="$1"
  if ! has_cmd python3; then
    return
  fi
  if [[ ! -f "${rust_file}" ]]; then
    return
  fi

  RUST_FILE="${rust_file}" python3 - <<'PY'
import json
import os

rust_file = os.environ.get("RUST_FILE", "")
ignore_csv = os.environ.get("IGNORE_DIRS_CSV", "")
ignore = [s.strip().lower() for s in ignore_csv.split(",") if s.strip()]

if not rust_file or not os.path.exists(rust_file):
    raise SystemExit(0)

def has_ignored_segment(value):
    if not value:
        return False
    val = str(value).lower().replace("\\", "/")
    for token in ignore:
        if f"/{token}/" in val or val.endswith(f"/{token}"):
            return True
    return False

def should_ignore(dep):
    if not isinstance(dep, dict):
        return True
    name = str(dep.get("name", "")).lower()
    for token in ignore:
        if token and token in name:
            return True
    for key in ("license_file", "notice_file", "repository", "manifest_path"):
        if has_ignored_segment(dep.get(key)):
            return True
    return False

with open(rust_file) as f:
    try:
        deps = json.load(f)
    except json.JSONDecodeError:
        deps = []

if isinstance(deps, list):
    filtered = [dep for dep in deps if not should_ignore(dep)]
else:
    filtered = deps

with open(rust_file, "w") as f:
    json.dump(filtered, f, indent=2)
PY
}

prepare_rust_manifest() {
  local manifest="$1"
  local crate_dir
  local tmp_dir
  local manifest_path="${manifest}"
  crate_dir="$(dirname "${manifest}")"

  if ! has_cmd python3; then
    echo "|${manifest_path}"
    return
  fi

  tmp_dir="$(mktemp -d)"
  if ! MANIFEST_SRC="${manifest}" MANIFEST_DST="${tmp_dir}/Cargo.toml" ROOT_DIR="${ROOT_DIR}" IGNORE_DIRS_CSV="${IGNORE_DIRS_CSV}" python3 - <<'PY'
import os
import re

src = os.environ.get("MANIFEST_SRC", "")
dst = os.environ.get("MANIFEST_DST", "")
root = os.environ.get("ROOT_DIR", "")
ignore_csv = os.environ.get("IGNORE_DIRS_CSV", "")
ignore = [s.strip().lower() for s in ignore_csv.split(",") if s.strip()]

if not src or not dst or not root:
    raise SystemExit(1)

base_dir = os.path.dirname(src)

with open(src) as f:
    lines = f.readlines()

def read_version(dir_name):
    cargo_toml = os.path.join(root, dir_name, "Cargo.toml")
    if not os.path.exists(cargo_toml):
        return None
    with open(cargo_toml) as f:
        in_package = False
        for line in f:
            stripped = line.strip()
            if stripped.startswith("[package]"):
                in_package = True
                continue
            if in_package and stripped.startswith("[") and stripped.endswith("]"):
                break
            if in_package:
                match = re.match(r'version\s*=\s*"([^"]+)"', stripped)
                if match:
                    return match.group(1)
    return None

dir_versions = {}
dir_paths = {}
for token in ignore:
    if not token:
        continue
    abs_dir = os.path.join(root, token)
    if os.path.exists(abs_dir):
        dir_paths[token] = abs_dir
        version = read_version(token)
        if version:
            dir_versions[token] = version

def resolve_path(path_value):
    abs_path = os.path.abspath(os.path.join(base_dir, path_value))
    if os.path.exists(abs_path):
        return abs_path
    normalized = path_value.replace("\\", "/").rstrip("/")
    for token in ignore:
        if not token:
            continue
        if normalized == token or normalized.endswith(f"/{token}"):
            candidate = os.path.join(root, token)
            if os.path.exists(candidate):
                return candidate
    return abs_path

out_lines = []
path_re = re.compile(r'(path\s*=\s*")([^"]+)(")')
table_re = re.compile(r'^\[([^\]]+)\]')
dep_line_re = re.compile(r'^\s*([A-Za-z0-9_.-]+)\s*=')
current_dep = None
skip_table = False
for line in lines:
    table_match = table_re.match(line.strip())
    if table_match:
        table_name = table_match.group(1)
        current_dep = None
        skip_table = False
        if "." in table_name:
            current_dep = table_name.rsplit(".", 1)[-1].strip()
            if current_dep in ignore:
                skip_table = True
                continue

    if skip_table:
        continue

    dep_match = dep_line_re.match(line)
    if dep_match:
        dep_name = dep_match.group(1).strip()
        if dep_name in ignore:
            continue

    match = path_re.search(line)
    if not match:
        if current_dep and current_dep in dir_versions:
            line = re.sub(
                r'version\s*=\s*"[^"]+"',
                f'version = "{dir_versions[current_dep]}"',
                line,
            )
        out_lines.append(line)
        continue
    resolved = resolve_path(match.group(2))
    updated_line = line[:match.start(2)] + resolved + line[match.end(2):]
    for token, abs_dir in dir_paths.items():
        if resolved == abs_dir or resolved.startswith(abs_dir + os.sep):
            if token in dir_versions:
                updated_line = re.sub(
                    r'version\s*=\s*"[^"]+"',
                    f'version = "{dir_versions[token]}"',
                    updated_line,
                )
            break
    out_lines.append(updated_line)

with open(dst, "w") as f:
    f.writelines(out_lines)
PY
  then
    rm -rf "${tmp_dir}"
    echo "|${manifest_path}"
    return
  fi

  if [[ -f "${crate_dir}/Cargo.lock" ]]; then
    cp "${crate_dir}/Cargo.lock" "${tmp_dir}/Cargo.lock"
  fi

  echo "${tmp_dir}|${tmp_dir}/Cargo.toml"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --install)
      INSTALL_MISSING=true
      shift
      ;;
    --scan)
      SCAN_TARGETS="${2:-}"
      shift 2
      ;;
    --fetch)
      FETCH_LICENSES=true
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

has_cmd() {
  command -v "$1" >/dev/null 2>&1
}

scan_enabled() {
  local target="$1"
  if [[ "${SCAN_TARGETS}" == "all" ]]; then
    return 0
  fi
  IFS=',' read -r -a targets <<< "${SCAN_TARGETS}"
  for t in "${targets[@]}"; do
    if [[ "${t}" == "${target}" ]]; then
      return 0
    fi
  done
  return 1
}

report_go() {
  local -a go_mods=()
  while IFS= read -r mod; do
    go_mods+=("${mod}")
  done < <(find_manifests "go.mod")
  if (( ${#go_mods[@]} == 0 )); then
    return
  fi

  if ! scan_enabled "go"; then
    return
  fi

  if ! has_cmd go; then
    echo "go not found; skipping Go license scan."
    return
  fi

  if has_cmd go-licenses; then
    for mod in "${go_mods[@]}"; do
      local mod_dir
      local slug
      mod_dir="$(dirname "${mod}")"
      slug="$(slugify_path "${mod_dir}")"
      echo "Scanning Go dependencies with go-licenses in ${mod_dir}..."
      (cd "${mod_dir}" && go-licenses csv ./... > "${OUT_DIR}/go-licenses-${slug}.csv")
      echo "Go licenses written to ${OUT_DIR}/go-licenses-${slug}.csv"
    done
  else
    if [[ "${INSTALL_MISSING}" == "true" ]]; then
      echo "Installing go-licenses..."
      (cd "${ROOT_DIR}" && go install github.com/google/go-licenses@latest)
    else
      echo "go-licenses not found. Install with:"
      echo "  go install github.com/google/go-licenses@latest"
      return
    fi
    if has_cmd go-licenses; then
      for mod in "${go_mods[@]}"; do
        local mod_dir
        local slug
        mod_dir="$(dirname "${mod}")"
        slug="$(slugify_path "${mod_dir}")"
        echo "Scanning Go dependencies with go-licenses in ${mod_dir}..."
        (cd "${mod_dir}" && go-licenses csv ./... > "${OUT_DIR}/go-licenses-${slug}.csv")
        echo "Go licenses written to ${OUT_DIR}/go-licenses-${slug}.csv"
      done
    else
      echo "go-licenses still not found after install attempt."
    fi
  fi
}

report_rust() {
  local -a cargo_manifests=()
  while IFS= read -r manifest; do
    cargo_manifests+=("${manifest}")
  done < <(find_manifests "Cargo.toml")
  if (( ${#cargo_manifests[@]} == 0 )); then
    return
  fi

  if ! scan_enabled "rust"; then
    return
  fi

  if ! has_cmd cargo; then
    echo "cargo not found; skipping Rust license scan."
    return
  fi

  if has_cmd cargo-license; then
    for manifest in "${cargo_manifests[@]}"; do
      local crate_dir
      local slug
      local temp_dir
      local manifest_path
      crate_dir="$(dirname "${manifest}")"
      slug="$(slugify_path "${crate_dir}")"
      echo "Scanning Rust dependencies with cargo-license in ${crate_dir}..."
      IFS='|' read -r temp_dir manifest_path < <(prepare_rust_manifest "${manifest}")
      (cd "${crate_dir}" && cargo license --current-dir "${crate_dir}" --manifest-path "${manifest_path}" --json > "${OUT_DIR}/rust-licenses-${slug}.json")
      if [[ -n "${temp_dir}" ]]; then
        rm -rf "${temp_dir}"
      fi
      filter_rust_licenses "${OUT_DIR}/rust-licenses-${slug}.json"
      echo "Rust licenses written to ${OUT_DIR}/rust-licenses-${slug}.json"
    done
  else
    if [[ "${INSTALL_MISSING}" == "true" ]]; then
      echo "Installing cargo-license..."
      (cd "${ROOT_DIR}" && cargo install cargo-license)
    else
      echo "cargo-license not found. Install with:"
      echo "  cargo install cargo-license"
      return
    fi
    if has_cmd cargo-license; then
      for manifest in "${cargo_manifests[@]}"; do
        local crate_dir
        local slug
        local temp_dir
        local manifest_path
        crate_dir="$(dirname "${manifest}")"
        slug="$(slugify_path "${crate_dir}")"
        echo "Scanning Rust dependencies with cargo-license in ${crate_dir}..."
        IFS='|' read -r temp_dir manifest_path < <(prepare_rust_manifest "${manifest}")
        (cd "${crate_dir}" && cargo license --current-dir "${crate_dir}" --manifest-path "${manifest_path}" --json > "${OUT_DIR}/rust-licenses-${slug}.json")
        if [[ -n "${temp_dir}" ]]; then
          rm -rf "${temp_dir}"
        fi
        filter_rust_licenses "${OUT_DIR}/rust-licenses-${slug}.json"
        echo "Rust licenses written to ${OUT_DIR}/rust-licenses-${slug}.json"
      done
    else
      echo "cargo-license still not found after install attempt."
    fi
  fi
}

report_go
report_rust

combine_reports() {
  local combined="${OUT_DIR}/all-licenses.csv"

  if ! has_cmd python3; then
    echo "python3 not found; skipping combined report."
    return
  fi

  python3 - <<'PY'
import csv
import glob
import os
import json

out_dir = os.environ.get("OUT_DIR", "")
combined = os.path.join(out_dir, "all-licenses.csv")
ignore_csv = os.environ.get("IGNORE_DIRS_CSV", "")
ignore = [s.strip().lower() for s in ignore_csv.split(",") if s.strip()]

def has_ignored_token(value):
    if not value:
        return False
    val = str(value).lower().replace("\\", "/")
    for token in ignore:
        if token in val:
            return True
    return False

def extract_license_value(row):
    for key in row.keys():
        key_l = key.lower()
        if key_l in ("license", "licenses"):
            return row[key]
    for key in row.keys():
        if "license" in key.lower():
            return row[key]
    return None

licenses = {}

for go_file in sorted(glob.glob(os.path.join(out_dir, "go-licenses*.csv"))):
    with open(go_file, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 3:
                continue
            module = (row[0] or "").strip()
            if has_ignored_token(module):
                continue
            lic = (row[2] or "").strip()
            if not lic or has_ignored_token(lic):
                continue
            licenses.setdefault(lic, set()).add("go")

for rust_file in sorted(glob.glob(os.path.join(out_dir, "rust-licenses*.json"))):
    with open(rust_file) as f:
        try:
            rust_deps = json.load(f)
        except json.JSONDecodeError:
            rust_deps = []
        if isinstance(rust_deps, list):
            for dep in rust_deps:
                if not isinstance(dep, dict):
                    continue
                name = dep.get("name") or ""
                if has_ignored_token(name):
                    continue
                lic = dep.get("license") or dep.get("licenses")
                if not lic:
                    continue
                lic = str(lic).strip()
                if not lic or has_ignored_token(lic):
                    continue
                licenses.setdefault(lic, set()).add("rust")

with open(combined, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["license", "languages"])
    for lic in sorted(licenses.keys()):
        langs = ",".join(sorted(licenses[lic]))
        writer.writerow([lic, langs])

print(f"Combined license report written to {combined}")
PY
}

fetch_go_licenses() {
  local cache_dir="${OUT_DIR}/cache/go"
  local overrides_file="${OUT_DIR}/license-overrides.csv"

  if [[ "${FETCH_LICENSES}" != "true" ]]; then
    return
  fi

  if ! scan_enabled "go"; then
    return
  fi

  if ! has_cmd curl; then
    echo "curl not found; skipping license fetch."
    return
  fi

  mkdir -p "${cache_dir}"

  python3 - <<'PY'
import csv
import glob
import os
import re
import subprocess
from urllib.parse import urlparse

out_dir = os.environ.get("OUT_DIR", "")
cache_dir = os.path.join(out_dir, "cache", "go")
overrides_file = os.path.join(out_dir, "license-overrides.csv")

def safe_name(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)

overrides = {}
if os.path.exists(overrides_file):
    with open(overrides_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            orig = (row.get("original_url") or "").strip()
            over = (row.get("override_url") or "").strip()
            if orig and over:
                overrides[orig] = over

for go_file in sorted(glob.glob(os.path.join(out_dir, "go-licenses*.csv"))):
    with open(go_file, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 2:
                continue
            module, url = row[0].strip(), row[1].strip()
            if not url.startswith("http"):
                continue
            url = overrides.get(url, url)
            parsed = urlparse(url)
            host = parsed.netloc.replace(":", "_")
            path = parsed.path.strip("/").replace("/", "_")
            fname = safe_name(f"{module}__{host}__{path}.txt")
            dest = os.path.join(cache_dir, fname)
            if os.path.exists(dest) and os.path.getsize(dest) > 0:
                continue
            try:
                subprocess.run(
                    ["curl", "-fsSL", url, "-o", dest],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                print(f"Fetched {url}")
            except subprocess.CalledProcessError:
                if os.path.exists(dest):
                    os.remove(dest)
                print(f"Failed to fetch {url}")
PY
}

fetch_go_licenses

combine_reports
