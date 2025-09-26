#!/usr/bin/env bash
set -euo pipefail

cd cli

# Enforce formatting
cargo fmt --all

# Lint everything (lib, bins, tests, benches, examples), treat warnings as errors
cargo clippy --fix --allow-dirty --all-targets --all-features --no-deps -- -D warnings
