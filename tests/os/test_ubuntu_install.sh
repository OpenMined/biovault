#!/bin/bash

set -e

echo "========================================="
echo "BioVault Ubuntu Installation Test"
echo "========================================="
echo ""

echo "Environment Info:"
echo "OS: $(lsb_release -ds 2>/dev/null || cat /etc/os-release | grep PRETTY_NAME | cut -d'=' -f2 | tr -d '\"')"
echo "Architecture: $(uname -m)"
echo ""

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to uninstall packages that might interfere
cleanup_existing() {
    echo "Cleaning up potentially conflicting packages..."

    # Remove existing Java installations if present (GitHub runners have Java pre-installed)
    # We'll remove it to test our installation process
    if command_exists java; then
        JAVA_VERSION=$(java -version 2>&1 | head -1 | grep -oP '\d+' | head -1)
        echo "Found existing Java installation (version $JAVA_VERSION), removing to test fresh install..."
        sudo apt-get remove -y --purge openjdk-* default-jdk* || true
        sudo apt-get autoremove -y || true
    fi

    # Remove existing Nextflow if present
    if [ -f /usr/local/bin/nextflow ]; then
        echo "Found existing Nextflow, removing..."
        sudo rm -f /usr/local/bin/nextflow
    fi

    # Docker is usually not installed on GitHub runners for Ubuntu, but check anyway
    if command_exists docker; then
        echo "Docker already installed, skipping removal to preserve runner state"
    fi

    echo ""
}

# Only cleanup on CI to avoid affecting local development
if [ -n "$CI" ] || [ -n "$GITHUB_ACTIONS" ]; then
    cleanup_existing
fi

echo "========================================="
echo "Testing bv check (before setup)"
echo "========================================="
./cli/target/release/bv check || true
echo ""

echo "========================================="
echo "Testing bv setup"
echo "========================================="
./cli/target/release/bv setup
echo ""

echo "========================================="
echo "Verifying installations"
echo "========================================="

# Check Java installation and verify it meets minimum version
if command_exists java; then
    JAVA_VERSION_OUTPUT=$(java -version 2>&1 | head -1)
    JAVA_VERSION=$(echo "$JAVA_VERSION_OUTPUT" | grep -oP '\d+' | head -1)
    echo "✓ Java installed: $JAVA_VERSION_OUTPUT"

    if [ "$JAVA_VERSION" -ge 17 ]; then
        echo "  Version $JAVA_VERSION meets minimum requirement (>=17)"
    else
        echo "✗ Java version $JAVA_VERSION is below minimum requirement (17)"
        exit 1
    fi
else
    echo "✗ Java not found"
    exit 1
fi

# Check Nextflow installation
if [ -f /usr/local/bin/nextflow ] || command_exists nextflow; then
    echo "✓ Nextflow installed"
    nextflow -version || /usr/local/bin/nextflow -version || true
else
    echo "✗ Nextflow not found"
    exit 1
fi

# Check Docker installation (might not be available on all CI runners)
if command_exists docker; then
    echo "✓ Docker installed:"
    docker --version
else
    echo "⚠️  Docker not installed (may not be available in CI environment)"
fi

# Check SyftBox installation
if command_exists syftbox; then
    echo "✓ SyftBox installed:"
    syftbox -v || true
else
    echo "⚠️  SyftBox not installed (setup-only mode)"
fi

echo ""
echo "========================================="
echo "Testing bv check (after setup)"
echo "========================================="
./cli/target/release/bv check
echo ""

echo "========================================="
echo "✓ Ubuntu installation test completed"
echo "========================================="