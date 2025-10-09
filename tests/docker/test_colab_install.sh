#!/bin/bash

# Script to test BioVault setup in a Colab-like environment
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

echo "========================================="
echo "BioVault Colab Installation Test"
echo "========================================="
echo ""

# Use pre-built image if specified (e.g., in CI), otherwise build locally
if [ -n "$DOCKER_IMAGE_NAME" ]; then
    echo "Using pre-built Docker image: $DOCKER_IMAGE_NAME"
    IMAGE_NAME="$DOCKER_IMAGE_NAME"
else
    echo "Building Docker image..."
    docker build -t biovault-colab-test "${SCRIPT_DIR}"
    IMAGE_NAME="biovault-colab-test"
fi

echo ""
echo "Running tests in container..."
echo ""

# Detect if running in CI
if [ -n "$CI" ] || [ -n "$GITHUB_ACTIONS" ]; then
    DOCKER_FLAGS="--rm"
    echo "Running in CI mode (non-interactive)"
else
    DOCKER_FLAGS="--rm -it"
    echo "Running in interactive mode"
fi

# Run the container with the source mounted
docker run $DOCKER_FLAGS \
    -v "${PROJECT_ROOT}/cli:/workspace/cli_src:ro" \
    "$IMAGE_NAME" \
    -c '
set -e

echo "========================================="
echo "Container Environment Info"
echo "========================================="
echo "OS Version: $(lsb_release -a 2>/dev/null | grep Description | cut -f2)"
echo "Java Version: $(java -version 2>&1 | head -1)"
echo "Rust Version: $(rustc --version)"
echo "COLAB_RELEASE_TAG: $COLAB_RELEASE_TAG"
echo "UV Version: $(uv --version 2>&1 || echo "uv not found")"
echo "PATH: $PATH"
echo "UV location: $(which uv 2>/dev/null || echo "uv not in PATH")"
echo ""

echo "========================================="
echo "Copying source and Building BioVault CLI"
echo "========================================="
cp -r /workspace/cli_src /workspace/cli
cd /workspace/cli
cargo build --release
echo "✓ Build successful"
echo ""

echo "========================================="
echo "Testing bv check (before setup)"
echo "========================================="
./target/release/bv check || true
echo ""

echo "========================================="
echo "Testing bv setup"
echo "========================================="
./target/release/bv setup
echo ""

echo "========================================="
echo "Verifying installations"
echo "========================================="

# Check Java upgrade
echo -n "Java version after setup: "
java -version 2>&1 | head -1

# Check Nextflow installation
if [ -f /usr/local/bin/nextflow ]; then
    echo "✓ Nextflow installed at /usr/local/bin/nextflow"
    /usr/local/bin/nextflow -version || true
else
    echo "✗ Nextflow not found at /usr/local/bin/nextflow"
fi

echo ""
echo "========================================="
echo "Testing bv check (after setup)"
echo "========================================="
./target/release/bv check
echo ""

echo "========================================="
echo "✓ All tests completed"
echo "========================================="
'