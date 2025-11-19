#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./decrypt.sh <relative-path-under-datasites> [--vault PATH] [--data-root DIR] [--shadow-root DIR] [--identity EMAIL] [--output FILE]
#
# Examples:
#   ./decrypt.sh client2@sandbox.local/app_data/biovault/rpc/message/abc.response
#   ./decrypt.sh client1@sandbox.local/shared/biovault/submissions/foo/bar.response --vault sandbox/client1@sandbox.local/.syc
#
# Defaults assume sandbox layout under the repo root:
#   vault:       sandbox/<email>/.syc       (derived from the first path component)
#   data-root:   sandbox/<email>/datasites  (derived from the first path component)
#   shadow-root: sandbox/<email>/unencrypted
#   identity:    <email> (derived from the first path component)
#   output:      stdout (omit --output to print)

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SYC_BIN=syc

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <relative-path-under-datasites> [--vault PATH] [--data-root DIR] [--shadow-root DIR] [--identity EMAIL] [--output FILE]" >&2
  exit 1
fi

INPUT_RAW="$1"; shift

# Normalize path and derive host (owner) email from the sandbox layout
trimmed="${INPUT_RAW#./}"
if [[ "$trimmed" == sandbox/* ]]; then
  trimmed="${trimmed#sandbox/}"
elif [[ "$trimmed" == */sandbox/* ]]; then
  trimmed="${trimmed#*sandbox/}"
fi

host="${trimmed%%/*}"
rel="$trimmed"
if [[ "$trimmed" == "$host/datasites/"* ]]; then
  rel="${trimmed#${host}/datasites/}"
elif [[ "$trimmed" == "$host/"* ]]; then
  rel="${trimmed#${host}/}"
fi

vault="sandbox/${host}/.syc"
data_root="sandbox/${host}/datasites"
shadow_root="sandbox/${host}/unencrypted"
identity="$host"
output=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vault) vault="$2"; shift ;;
    --data-root) data_root="$2"; shift ;;
    --shadow-root) shadow_root="$2"; shift ;;
    --identity) identity="$2"; shift ;;
    --output) output="$2"; shift ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
  shift
done

cmd=("$SYC_BIN" bytes read
  --vault "$vault"
  --data-root "$data_root"
  --shadow-root "$shadow_root"
  --identity "$identity"
  --relative "$rel"
)

emit_to_file=false
if [[ -n "$output" ]]; then
  cmd+=( --output "$output" )
  emit_to_file=true
fi

echo "Running: ${cmd[*]}" >&2

if $emit_to_file; then
  "${cmd[@]}"
  content="$(cat "$output")"
else
  content="$("${cmd[@]}")"
  printf '%s\n' "$content"
fi

# Decode JSON body if present (best-effort). Strip non-JSON lines (e.g., syc logs) first.
python3 <<PY
import base64, json
content = """$content"""
start = content.find("{")
if start == -1:
    import sys
    sys.exit(0)
json_text = content[start:]
try:
    obj = json.loads(json_text)
except Exception:
    import sys
    sys.exit(0)

body = obj.get("body")
if body is None:
    import sys
    sys.exit(0)

try:
    raw = body.encode() if isinstance(body, str) else body
    decoded = base64.b64decode(raw).decode("utf-8", errors="replace")
    print("\n--- Decoded body ---")
    print(decoded)
except Exception as e:
    print(f"\n--- Body (decode error: {e}) ---")
    print(body)
PY
