#!/usr/bin/env bash
# Quick test for Nextflow via Podman nested containers
# NO CLI needed - just tests that Podman can spawn containers from within containers
#
# Usage: ./tests/scripts/quick-nextflow-test.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BIOSCRIPT_DIR="$(cd "$ROOT_DIR/../bioscript" 2>/dev/null || echo "$ROOT_DIR/../bioscript")"

echo "=== Quick Nextflow + Podman Test ==="
echo "Testing that Nextflow can spawn task containers via Podman socket"
echo ""

# Create minimal test workflow
TEST_DIR="/tmp/quick-nf-test-$$"
mkdir -p "$TEST_DIR"
echo "Test directory: $TEST_DIR"

cat > "$TEST_DIR/test.nf" <<'EOF'
nextflow.enable.dsl=2

process HELLO {
    container 'ghcr.io/openmined/bioscript:latest'

    output:
        stdout

    script:
    """
    echo "Hello from nested container!"
    bioscript --help | head -5
    """
}

workflow {
    HELLO() | view
}
EOF

cat > "$TEST_DIR/nextflow.config" <<'EOF'
process.executor = 'local'
podman.enabled = true
process.shell = ['/bin/sh', '-ue']
EOF

# Get Podman socket path
PODMAN_SOCKET="${XDG_RUNTIME_DIR:-/run/user/$(id -u)}/podman/podman.sock"
if [[ ! -S "$PODMAN_SOCKET" ]]; then
    echo "ERROR: Podman socket not found at $PODMAN_SOCKET"
    echo "Is Podman machine running?"
    exit 1
fi

echo "Podman socket: $PODMAN_SOCKET"

# Convert paths for container mounts
# Git Bash /c/... -> /mnt/c/...
container_path() {
    echo "$1" | sed 's|^/\([a-zA-Z]\)/|/mnt/\1/|'
}

TEST_DIR_MOUNT=$(container_path "$TEST_DIR")

echo ""
echo "=== Running Nextflow in container ==="
podman run --rm \
    --userns=keep-id \
    --security-opt label=disable \
    -e "NXF_HOME=/tmp/.nextflow" \
    -v "$PODMAN_SOCKET:/run/podman/podman.sock" \
    -e "CONTAINER_HOST=unix:///run/podman/podman.sock" \
    -v "$TEST_DIR_MOUNT:$TEST_DIR_MOUNT" \
    -w "$TEST_DIR_MOUNT" \
    ghcr.io/openmined/nextflow-runner:25.10.2 \
    nextflow run test.nf -c nextflow.config

echo ""
echo "=== Test Complete ==="
echo "Podman nested containers working!"

# Cleanup
rm -rf "$TEST_DIR"
