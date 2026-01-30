#!/usr/bin/env bash
set -euo pipefail

# GitHub Actions Windows runners often provide `python` but not `python3` on PATH.
# Normalize so the rest of the script can keep using `python3`.
if ! command -v python3 >/dev/null 2>&1; then
  if command -v python >/dev/null 2>&1; then
    python3() { python "$@"; }
  else
    echo "Missing required tool: python3 (or python)" >&2
    exit 1
  fi
fi

usage() {
  cat <<'EOF' >&2
Usage: ./test-scenario.sh [options] <scenario.yaml>

Options:
  --client-mode MODE   SyftBox client: go|rust|mixed|embedded (default: rust)
  --sandbox DIR        Sandbox root (default: ./sandbox)
  --rust-client-bin P  Path to Rust client binary (optional)
  --skip-rust-build    Do not build Rust client (requires binary exists)
  --no-reset           Do not reset devstack/sandbox (reuse existing state)
  --force              Force scenario steps to re-run (overrides skip-done)
  --allele-count N     Override allele-freq synthetic file count (default: 10)
  --he                 Use HE aggregation path in syqure-flow (manual ciphertext exchange)
  --docker             Force Docker mode for syqure runtime
  --podman             Force Podman runtime (sets BIOVAULT_CONTAINER_RUNTIME=podman)
  --keep-containers    Keep syqure containers on failure (for logs/debugging)
  -h, --help           Show this message

Examples:
  ./test-scenario.sh tests/scenarios/inbox-ping-pong.yaml
  ./test-scenario.sh --client-mode go tests/scenarios/inbox-ping-pong.yaml
  ./test-scenario.sh --sandbox sandbox-rs tests/scenarios/inbox-ping-pong.yaml
  ./test-scenario.sh --client-mode embedded tests/scenarios/inbox-ping-pong.yaml
  ./test-scenario.sh --docker tests/scenarios/syqure-distributed.yaml
  ./test-scenario.sh --podman tests/scenarios/syqure-distributed.yaml
  ./test-scenario.sh --podman --keep-containers tests/scenarios/syqure-distributed.yaml
  ./test-scenario.sh --no-reset tests/scenarios/allele-freq-syqure.yaml
EOF
}

CLIENT_MODE="rust"
SANDBOX_DIR=""
RUST_CLIENT_BIN=""
SKIP_RUST_BUILD=0
USE_DOCKER=0
USE_PODMAN=0
KEEP_CONTAINERS=0
NO_RESET=0
FORCE_RUN=0
ALLELE_COUNT=""
USE_HE=0
SCENARIO=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --client-mode)
      [[ $# -lt 2 ]] && { echo "Missing value for --client-mode" >&2; usage; exit 1; }
      CLIENT_MODE="${2:-}"
      shift
      ;;
    --sandbox)
      [[ $# -lt 2 ]] && { echo "Missing value for --sandbox" >&2; usage; exit 1; }
      SANDBOX_DIR="${2:-}"
      shift
      ;;
    --rust-client-bin)
      [[ $# -lt 2 ]] && { echo "Missing value for --rust-client-bin" >&2; usage; exit 1; }
      RUST_CLIENT_BIN="${2:-}"
      shift
      ;;
    --skip-rust-build)
      SKIP_RUST_BUILD=1
      ;;
    --no-reset)
      NO_RESET=1
      ;;
    --force)
      FORCE_RUN=1
      ;;
    --allele-count)
      [[ $# -lt 2 ]] && { echo "Missing value for --allele-count" >&2; usage; exit 1; }
      ALLELE_COUNT="${2:-}"
      shift
      ;;
    --he)
      USE_HE=1
      ;;
    --docker)
      USE_DOCKER=1
      ;;
    --podman)
      USE_PODMAN=1
      USE_DOCKER=1
      ;;
    --keep-containers)
      KEEP_CONTAINERS=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      SCENARIO="$1"
      ;;
  esac
  shift
done

if [[ -z "$SCENARIO" ]]; then
  usage
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if (( USE_HE )); then
  if [[ -f "$SCENARIO" ]]; then
    if grep -q "syqure-flow/flow.yaml" "$SCENARIO"; then
      mkdir -p "$ROOT_DIR/logs"
      SCENARIO_HE="$ROOT_DIR/logs/$(basename "$SCENARIO" .yaml).he.yaml"
      sed 's#syqure-flow/flow.yaml#syqure-flow/flow-he.yaml#g' "$SCENARIO" >"$SCENARIO_HE"
      SCENARIO="$SCENARIO_HE"
    else
      echo "Warning: --he set but scenario does not reference syqure-flow/flow.yaml" >&2
    fi
  fi
fi

echo "Building BioVault CLI (release)..."
(cd "$ROOT_DIR/cli" && cargo build --release)

echo "Running scenario: $SCENARIO"

# Preserve user JAVA_HOME/PATH for downstream shells (after sbenv activation)
export SCENARIO_JAVA_HOME="${JAVA_HOME:-}"
export SCENARIO_USER_PATH="$PATH"
if command -v nextflow >/dev/null 2>&1; then
  export BIOVAULT_BUNDLED_NEXTFLOW="$(command -v nextflow)"
fi

if [[ -n "$SANDBOX_DIR" ]]; then
  export SANDBOX_DIR="$SANDBOX_DIR"
fi

if (( USE_PODMAN )); then
  export BIOVAULT_CONTAINER_RUNTIME="podman"
  export CONTAINERS_MACHINE_PROVIDER="${CONTAINERS_MACHINE_PROVIDER:-hyperv}"
  export BIOVAULT_HYPERV_MOUNT="${BIOVAULT_HYPERV_MOUNT:-1}"
fi
if (( KEEP_CONTAINERS )); then
  export BIOVAULT_SYQURE_KEEP_CONTAINERS="1"
fi
if (( USE_DOCKER )); then
  # Ensure syqure.yaml switches into Docker mode when --docker/--podman is used.
  export SEQURE_MODE="${SEQURE_MODE:-docker}"
fi

if (( NO_RESET )); then
  export BV_DEVSTACK_RESET=0
fi
if (( FORCE_RUN )); then
  export BV_FLOW_FORCE=1
fi
if [[ -n "$ALLELE_COUNT" ]]; then
  export ALLELE_FREQ_COUNT="$ALLELE_COUNT"
fi
if (( USE_HE )); then
  export BV_SYQURE_HE=1
fi

# Let tests/scripts/devstack.sh decide which syftbox client to run.
export BV_DEVSTACK_CLIENT_MODE="$CLIENT_MODE"
if [[ -n "$RUST_CLIENT_BIN" ]]; then
  export BV_DEVSTACK_RUST_CLIENT_BIN="$RUST_CLIENT_BIN"
fi
if (( SKIP_RUST_BUILD )); then
  export BV_DEVSTACK_SKIP_RUST_BUILD=1
fi

if [[ "${CLIENT_MODE}" == "embedded" ]]; then
  # Ensure BioVault-hosted SyftBox runs in embedded mode when started by devstack.sh.
  export BV_SYFTBOX_BACKEND=embedded
fi

# Syqure runtime setup - only if scenario needs it (and not on Windows)
# Detect Windows (Git Bash / MSYS / Cygwin)
IS_WINDOWS=0
if [[ "$(uname -s)" == MINGW* ]] || [[ "$(uname -s)" == MSYS* ]] || [[ "$(uname -s)" == CYGWIN* ]]; then
  IS_WINDOWS=1
fi

# Check if scenario needs syqure or container runtime
NEEDS_SYQURE=0
if grep -qiE '(syqure|mpc|sequre)' "$SCENARIO" 2>/dev/null; then
  NEEDS_SYQURE=1
fi

NEEDS_CONTAINER=0
if grep -qiE '(syqure|mpc|sequre|nextflow|workflow\.nf|bv submit)' "$SCENARIO" 2>/dev/null; then
  NEEDS_CONTAINER=1
fi

if (( IS_WINDOWS )) && (( NEEDS_CONTAINER )); then
  RUNTIME_PREF="${BIOVAULT_CONTAINER_RUNTIME:-}"
  if (( USE_PODMAN )); then
    RUNTIME_PREF="podman"
  elif (( USE_DOCKER )); then
    RUNTIME_PREF="docker"
  fi

  check_docker() {
    command -v docker >/dev/null 2>&1 || return 1
    docker info >/dev/null 2>&1
  }

  check_podman() {
    command -v podman >/dev/null 2>&1 || return 1
    podman info >/dev/null 2>&1
  }

  if [[ "$RUNTIME_PREF" == "podman" ]] || { [[ -z "$RUNTIME_PREF" ]] && check_podman && ! check_docker; }; then
    export BIOVAULT_CONTAINER_RUNTIME="podman"
    export CONTAINERS_MACHINE_PROVIDER="${CONTAINERS_MACHINE_PROVIDER:-hyperv}"
    # Hyper-V mounts are unreliable on Windows; force VM-copy mode.
    export BIOVAULT_HYPERV_MOUNT="0"
  fi

  case "$RUNTIME_PREF" in
    podman)
      if ! check_podman; then
        echo "Podman is required for this scenario but is not running. Start it with 'podman machine start'." >&2
        exit 1
      fi
      ;;
    docker)
      if ! check_docker; then
        echo "Docker is required for this scenario but is not running. Start Docker Desktop." >&2
        exit 1
      fi
      ;;
    *)
      if check_docker; then
        : # ok
      elif check_podman; then
        : # ok
      else
        echo "This scenario requires a container runtime (Docker or Podman), but neither is running." >&2
        exit 1
      fi
      ;;
  esac
fi

if (( NEEDS_SYQURE )); then
  # Determine syqure directory
  if [[ -n "${BV_SYQURE_DIR:-}" ]]; then
    SYQURE_DIR="$BV_SYQURE_DIR"
  elif [[ -d "$ROOT_DIR/../syqure" ]]; then
    SYQURE_DIR="$ROOT_DIR/../syqure"
  else
    SYQURE_DIR="$ROOT_DIR/syqure"
  fi
  SYQURE_BIN_DEBUG="$SYQURE_DIR/target/debug/syqure"
  SYQURE_BIN_RELEASE="$SYQURE_DIR/target/release/syqure"
  SYQURE_BIN="$SYQURE_BIN_DEBUG"

  if (( IS_WINDOWS )); then
    # Syqure can't build on Windows - use Docker mode
    export BV_SYQURE_USE_DOCKER=1
    echo "Syqure mode: Docker (Windows - native build not supported)"
  elif (( USE_DOCKER )); then
    export BV_SYQURE_USE_DOCKER=1
    echo "Syqure mode: Docker"
  else
    # Preflight: if no bundle is available for native syqure, fall back to Docker.
    BUNDLE_OK=0
    HAS_PREBUILT_CODON=0
    if [[ -n "${SYQURE_BUNDLE_FILE:-}" && -f "${SYQURE_BUNDLE_FILE}" ]]; then
      BUNDLE_OK=1
    else
      if command -v rustc >/dev/null 2>&1; then
        HOST_TRIPLE="$(rustc -vV | sed -n 's/^host: //p')"
        if [[ -n "$HOST_TRIPLE" && -f "$SYQURE_DIR/bundles/${HOST_TRIPLE}.tar.zst" ]]; then
          BUNDLE_OK=1
        fi
        if [[ -n "$HOST_TRIPLE" && -f "$SYQURE_DIR/syqure/bundles/${HOST_TRIPLE}.tar.zst" ]]; then
          BUNDLE_OK=1
        fi
      fi
      if [[ -d "$SYQURE_DIR/bin/macos-arm64/codon" || -d "$SYQURE_DIR/bin/macos-x86_64/codon" || -d "$SYQURE_DIR/bin/linux-x86/codon" || -d "$SYQURE_DIR/bin/linux-arm64/codon" ]]; then
        HAS_PREBUILT_CODON=1
      fi
      if [[ -d "$ROOT_DIR/../codon/install/lib/codon" ]]; then
        BUNDLE_OK=1
      fi
    fi

    if [[ "${SYQURE_FORCE_BUNDLE_REBUILD:-0}" == "1" ]]; then
      BUNDLE_OK=0
    fi

    if (( ! BUNDLE_OK )); then
      if [[ -d "$SYQURE_DIR/.git" && -f "$SYQURE_DIR/.gitmodules" ]]; then
        echo "Syqure bundle/assets missing; initializing submodules..."
        (cd "$SYQURE_DIR" && git submodule update --init --recursive)
      fi

      if (( HAS_PREBUILT_CODON )) && [[ -x "$SYQURE_DIR/syqure_bins.sh" ]]; then
        echo "Building syqure bundle from prebuilts (syqure_bins.sh)..."
        (cd "$SYQURE_DIR" && ./syqure_bins.sh)
      elif [[ -x "$SYQURE_DIR/build_libs.sh" ]]; then
        echo "Building syqure bundle (build_libs.sh)..."
        (cd "$SYQURE_DIR" && ./build_libs.sh)
      fi

      # Re-check after attempted build.
      if [[ -n "${SYQURE_BUNDLE_FILE:-}" && -f "${SYQURE_BUNDLE_FILE}" ]]; then
        BUNDLE_OK=1
      else
        if command -v rustc >/dev/null 2>&1; then
          HOST_TRIPLE="$(rustc -vV | sed -n 's/^host: //p')"
          if [[ -n "$HOST_TRIPLE" && -f "$SYQURE_DIR/bundles/${HOST_TRIPLE}.tar.zst" ]]; then
            BUNDLE_OK=1
          fi
          if [[ -n "$HOST_TRIPLE" && -f "$SYQURE_DIR/syqure/bundles/${HOST_TRIPLE}.tar.zst" ]]; then
            BUNDLE_OK=1
          fi
        fi
      fi

      if (( ! BUNDLE_OK )); then
        USE_DOCKER=1
        export BV_SYQURE_USE_DOCKER=1
        echo "Syqure bundle not found; falling back to Docker. Set SYQURE_BUNDLE_FILE to use native."
      fi
    fi
  fi

  if (( USE_DOCKER )); then
    export BV_SYQURE_USE_DOCKER=1
    echo "Syqure mode: Docker"
  else
    # Native mode - build syqure if needed
    if [[ ! -x "$SYQURE_BIN" ]]; then
      if [[ -d "$SYQURE_DIR" ]]; then
        echo "Building syqure native binary (CI-style)..."
        # Mirror syqure CI smoke build: use precompiled Codon from bin/<platform>/codon if present.
        SYQURE_PLATFORM=""
        OS_NAME="$(uname -s | tr '[:upper:]' '[:lower:]')"
        ARCH_NAME="$(uname -m | tr '[:upper:]' '[:lower:]')"
        case "$OS_NAME" in
          darwin) OS_LABEL="macos" ;;
          linux) OS_LABEL="linux" ;;
          *) OS_LABEL="$OS_NAME" ;;
        esac
        case "$ARCH_NAME" in
          arm64|aarch64) ARCH_LABEL="arm64" ;;
          x86_64|amd64|i386|i686) ARCH_LABEL="x86" ;;
          *) ARCH_LABEL="$ARCH_NAME" ;;
        esac
        SYQURE_PLATFORM="${OS_LABEL}-${ARCH_LABEL}"
        BIN_ROOT="$SYQURE_DIR/bin/$SYQURE_PLATFORM/codon"
        if [[ -d "$BIN_ROOT" ]]; then
          export SYQURE_CPP_INCLUDE="$BIN_ROOT/include"
          export SYQURE_CPP_LIB_DIRS="$BIN_ROOT/lib/codon"
          # Fix broken libgmp symlink in precompiled bundle if needed (macOS).
          if [[ "$OS_LABEL" == "macos" ]]; then
            GMP_FILE="$BIN_ROOT/lib/codon/libgmp.dylib"
            if [[ -L "$GMP_FILE" && ! -e "$GMP_FILE" ]]; then
              echo "Repairing broken libgmp symlink in $BIN_ROOT..."
              GMP_SRC=""
              if command -v brew >/dev/null 2>&1; then
                GMP_PREFIX="$(brew --prefix gmp 2>/dev/null || true)"
                if [[ -n "$GMP_PREFIX" && -f "$GMP_PREFIX/lib/libgmp.dylib" ]]; then
                  GMP_SRC="$GMP_PREFIX/lib/libgmp.dylib"
                fi
              fi
              if [[ -z "$GMP_SRC" ]]; then
                for candidate in /opt/homebrew/opt/gmp/lib/libgmp.dylib /usr/local/opt/gmp/lib/libgmp.dylib; do
                  if [[ -f "$candidate" ]]; then
                    GMP_SRC="$candidate"
                    break
                  fi
                done
              fi
              if [[ -n "$GMP_SRC" ]]; then
                rm -f "$BIN_ROOT/lib/codon/libgmp.dylib" "$BIN_ROOT/lib/codon/libgmp.so"
                cp -L "$GMP_SRC" "$BIN_ROOT/lib/codon/libgmp.dylib"
                cp -L "$GMP_SRC" "$BIN_ROOT/lib/codon/libgmp.so"
              else
                echo "libgmp.dylib not found; install gmp or set SEQURE_GMP_PATH" >&2
                exit 1
              fi
            fi
          fi
        fi
        (cd "$SYQURE_DIR" && cargo build -p syqure) || {
          echo "Failed to build syqure. Use --docker flag for Docker mode." >&2
          exit 1
        }
      else
        echo "Syqure directory not found at $SYQURE_DIR. Use --docker flag for Docker mode." >&2
        exit 1
      fi
    fi
    # Prefer release binary for runtime stability unless explicitly disabled.
    if [[ "${BV_SYQURE_PREFER_RELEASE:-1}" == "1" ]]; then
      if [[ ! -x "$SYQURE_BIN_RELEASE" && -d "$SYQURE_DIR" ]]; then
        echo "Building syqure native binary (release)..."
        (cd "$SYQURE_DIR" && cargo build -p syqure --release) || {
          echo "Failed to build syqure (release). Falling back to debug binary." >&2
        }
      fi
      if [[ -x "$SYQURE_BIN_RELEASE" ]]; then
        export SEQURE_NATIVE_BIN="$SYQURE_BIN_RELEASE"
      fi
    fi

    if [[ -z "${SEQURE_NATIVE_BIN:-}" && -x "$SYQURE_BIN_DEBUG" ]]; then
      export SEQURE_NATIVE_BIN="$SYQURE_BIN_DEBUG"
    fi

    if [[ -x "${SEQURE_NATIVE_BIN:-}" ]]; then
      echo "Syqure mode: Native ($SEQURE_NATIVE_BIN)"
    fi
  fi
else
  echo "Syqure mode: Not needed for this scenario"
fi

# Syqure runs can take longer; increase background timeout unless user set it.
if (( NEEDS_SYQURE )); then
  if (( IS_WINDOWS )); then
    export SCENARIO_BG_TIMEOUT="${SCENARIO_BG_TIMEOUT:-3600}"
  else
    export SCENARIO_BG_TIMEOUT="${SCENARIO_BG_TIMEOUT:-1200}"
  fi
  echo "Background process timeout: ${SCENARIO_BG_TIMEOUT}s"
  if (( IS_WINDOWS )); then
    export SCENARIO_TAIL_SYQURE_LOGS="${SCENARIO_TAIL_SYQURE_LOGS:-1}"
  fi
fi

TAIL_PID=""
PROGRESS_PID=""
start_syqure_log_tail() {
  if [[ "${SCENARIO_TAIL_SYQURE_LOGS:-0}" != "1" ]]; then
    return
  fi
  if (( ! IS_WINDOWS )); then
    return
  fi
  echo "Tailing syqure file_transport.log (Windows)..."
  powershell.exe -NoProfile -Command '& { $ErrorActionPreference = "SilentlyContinue"; $root = (Get-Location).Path + "\\sandbox"; $pos = @{}; while ($true) { Get-ChildItem -Recurse -Filter file_transport.log $root -ErrorAction SilentlyContinue | ForEach-Object { $p = $_.FullName; if (-not $pos.ContainsKey($p)) { $pos[$p] = 0 }; try { $fs = [IO.File]::Open($p, [IO.FileMode]::Open, [IO.FileAccess]::Read, [IO.FileShare]::ReadWrite); $fs.Seek([long]$pos[$p], [IO.SeekOrigin]::Begin) | Out-Null; $sr = New-Object IO.StreamReader($fs); while (-not $sr.EndOfStream) { $line = $sr.ReadLine(); if ($line) { Write-Host ("[syqure-log] {0}: {1}" -f $p, $line) } }; $pos[$p] = $fs.Position; $sr.Close(); $fs.Close() } catch {} }; Start-Sleep -Seconds 2 } }' &
  TAIL_PID=$!
}

stop_syqure_log_tail() {
  if [[ -n "${TAIL_PID:-}" ]]; then
    kill "$TAIL_PID" >/dev/null 2>&1 || true
  fi
  if [[ -n "${PROGRESS_PID:-}" ]]; then
    kill "$PROGRESS_PID" >/dev/null 2>&1 || true
  fi
}

start_syqure_progress_tail() {
  if [[ "${SCENARIO_TAIL_SYQURE_LOGS:-0}" != "1" ]]; then
    return
  fi
  if (( ! IS_WINDOWS )); then
    return
  fi
  echo "Tracking syqure message counts (Windows)..."
  powershell.exe -NoProfile -Command '& { $ErrorActionPreference = "SilentlyContinue"; $root = (Get-Location).Path + "\\sandbox"; while ($true) { $now = Get-Date -Format "HH:mm:ss"; $groups = Get-ChildItem -Recurse -Filter *.request $root -ErrorAction SilentlyContinue | Group-Object DirectoryName | Sort-Object Count -Descending | Select-Object -First 5; if ($groups) { foreach ($g in $groups) { Write-Host ("[syqure-progress] {0} {1} files in {2}" -f $now, $g.Count, $g.Name) } } else { Write-Host ("[syqure-progress] {0} no request files yet" -f $now) }; Start-Sleep -Seconds 10 } }' &
  PROGRESS_PID=$!
}

if python3 -c 'import yaml' >/dev/null 2>&1; then
  start_syqure_log_tail
  start_syqure_progress_tail
  python3 "$ROOT_DIR/scripts/run_scenario.py" "$SCENARIO"
  status=$?
  stop_syqure_log_tail
  exit $status
fi

start_syqure_log_tail
start_syqure_progress_tail
python3 "$ROOT_DIR/scripts/run_scenario.py" "$SCENARIO"
status=$?
stop_syqure_log_tail
exit $status
