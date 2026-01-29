#!/usr/bin/env bash
set -euo pipefail
uv run "$(dirname "$0")/build_master.py"
