#!/bin/bash
# Quick test script for Nextflow pipeline with Podman nested containers
# Usage: ./tests/scripts/test-nextflow-podman.sh [workflow_dir]

set -e

WORKFLOW_DIR="${1:-/c/Users/admin/AppData/Local/Temp/bvtest}"

# Convert Windows path to Podman/WSL format (/mnt/c/...)
if [[ "$WORKFLOW_DIR" == /c/* ]]; then
    PODMAN_PATH="/mnt${WORKFLOW_DIR}"
elif [[ "$WORKFLOW_DIR" == C:* ]] || [[ "$WORKFLOW_DIR" == c:* ]]; then
    # Convert C:\path to /mnt/c/path
    PODMAN_PATH="/mnt/c${WORKFLOW_DIR:2}"
    PODMAN_PATH="${PODMAN_PATH//\\//}"
else
    PODMAN_PATH="$WORKFLOW_DIR"
fi

echo "=== Nextflow Podman Test ==="
echo "Workflow dir: $WORKFLOW_DIR"
echo "Podman path:  $PODMAN_PATH"
echo ""

# Find Podman socket
PODMAN_SOCKET=""
if [[ -S /run/user/1000/podman/podman.sock ]]; then
    PODMAN_SOCKET="/run/user/1000/podman/podman.sock"
elif [[ -S "$XDG_RUNTIME_DIR/podman/podman.sock" ]]; then
    PODMAN_SOCKET="$XDG_RUNTIME_DIR/podman/podman.sock"
else
    echo "Error: Could not find Podman socket"
    echo "Make sure Podman machine is running: podman machine start"
    exit 1
fi

echo "Podman socket: $PODMAN_SOCKET"
echo ""

# Clean up old .nextflow and work directories that might have wrong permissions
if [[ -d "$WORKFLOW_DIR/.nextflow" ]]; then
    echo "Cleaning old .nextflow directory..."
    rm -rf "$WORKFLOW_DIR/.nextflow" 2>/dev/null || true
fi
if [[ -d "$WORKFLOW_DIR/work" ]]; then
    echo "Cleaning old work directory..."
    rm -rf "$WORKFLOW_DIR/work" 2>/dev/null || true
fi

# Create fresh directories
mkdir -p "$WORKFLOW_DIR/.nextflow"
mkdir -p "$WORKFLOW_DIR/work"

echo "=== Running Nextflow with Podman ==="
echo ""

# Run Nextflow in container with Podman socket mounted for nested containers
podman run --rm \
    --userns=keep-id \
    --security-opt label=disable \
    -e "NXF_HOME=/tmp/.nextflow" \
    -v "${PODMAN_SOCKET}:/run/podman/podman.sock" \
    -e "CONTAINER_HOST=unix:///run/podman/podman.sock" \
    -v "${PODMAN_PATH}:${PODMAN_PATH}:rw" \
    -w "${PODMAN_PATH}" \
    nextflow/nextflow:25.10.2 \
    nextflow run workflow.nf \
        -with-podman \
        -work-dir "${PODMAN_PATH}/work"

echo ""
echo "=== Test Complete ==="
echo "Check output in: $WORKFLOW_DIR"
