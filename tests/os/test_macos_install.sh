#!/bin/bash

set -e

echo "========================================="
echo "BioVault macOS Installation Test"
echo "========================================="
echo ""

echo "Environment Info:"
echo "OS: $(sw_vers -productName) $(sw_vers -productVersion)"
echo "Architecture: $(uname -m)"
echo ""

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to uninstall packages that might interfere
cleanup_existing() {
    echo "Cleaning up potentially conflicting packages..."

    # On GitHub runners, Homebrew is pre-installed, but we might need to clean up some packages
    if command_exists brew; then
        echo "Homebrew is installed"

        # Remove existing Java installations if present via Homebrew
        # Check for any openjdk installations
        echo "Checking for existing OpenJDK installations..."
        for pkg in $(brew list --formula | grep openjdk); do
            echo "Found $pkg via Homebrew, removing..."
            brew uninstall --ignore-dependencies "$pkg" 2>/dev/null || true
        done

        # Remove existing Nextflow if installed via Homebrew
        if brew list --formula | grep -q nextflow; then
            echo "Found existing Nextflow via Homebrew, removing..."
            brew uninstall --ignore-dependencies nextflow || true
        fi
    fi

    # Remove manually installed Nextflow
    if [ -f /usr/local/bin/nextflow ]; then
        echo "Found existing Nextflow at /usr/local/bin, removing..."
        sudo rm -f /usr/local/bin/nextflow || rm -f /usr/local/bin/nextflow
    fi

    echo ""
}

# Only cleanup on CI to avoid affecting local development
if [ -n "$CI" ] || [ -n "$GITHUB_ACTIONS" ]; then
    cleanup_existing
fi

# Verify bv is installed
if ! command_exists bv; then
    echo "Error: bv command not found in PATH"
    echo "PATH: $PATH"
    exit 1
fi

echo "Using bv binary at: $(which bv)"
echo "bv version: $(bv --version || echo 'unknown')"
echo ""

echo "========================================="
echo "Testing bv check (before setup)"
echo "========================================="
bv check || true
echo ""

echo "========================================="
echo "Testing bv setup"
echo "========================================="
bv setup
echo ""

echo "========================================="
echo "Verifying installations"
echo "========================================="

# Check Java installation
if command_exists java; then
    echo "✓ Java installed:"
    java -version 2>&1 | head -1
else
    echo "✗ Java not found"
    exit 1
fi

# Check Nextflow installation
if command_exists nextflow || [ -f /usr/local/bin/nextflow ]; then
    echo "✓ Nextflow installed"
    nextflow -version || /usr/local/bin/nextflow -version || true
else
    echo "✗ Nextflow not found"
    exit 1
fi

# Check Docker installation (Docker Desktop might not be running on CI)
if command_exists docker; then
    echo "✓ Docker command available"
    docker --version || echo "Docker command exists but daemon might not be running"
else
    echo "⚠️  Docker not installed (Docker Desktop installation requires GUI interaction)"
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
bv check
echo ""

echo "========================================="
echo "✓ macOS installation test completed"
echo "========================================="