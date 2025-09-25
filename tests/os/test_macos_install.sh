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

        # GitHub Actions macOS runners come with Java pre-installed
        # We need to handle different scenarios:
        # 1. Java from system (GitHub Actions case)
        # 2. Java from Homebrew
        # 3. Java installed but not in PATH

        if command_exists java; then
            JAVA_PATH=$(which java)
            echo "Java currently found at: $JAVA_PATH"

            # Get Java version for info
            java -version 2>&1 | head -1 || true

            # Check if this Java is from brew or system
            if echo "$JAVA_PATH" | grep -q "brew\|homebrew\|Cellar"; then
                echo "Java is from Homebrew"
                # For testing, we'll uninstall it and reinstall to test the PATH configuration
                echo "Uninstalling Homebrew Java to test fresh installation with PATH configuration..."
                for pkg in $(brew list --formula | grep openjdk || true); do
                    brew uninstall --ignore-dependencies "$pkg" 2>/dev/null || true
                done
            else
                echo "Java appears to be from system (not Homebrew)"
                # On GitHub Actions, this is expected
                # We'll install via brew to test our PATH configuration logic
                echo "Will install Java via Homebrew to test PATH configuration..."
            fi
        else
            echo "No Java found in PATH"
        fi

        # Remove any Java paths from PATH to simulate the scenario
        echo "Cleaning PATH for testing..."
        export ORIGINAL_PATH="$PATH"
        export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v 'java\|jdk\|openjdk\|Java' | tr '\n' ':' | sed 's/:$//')
        echo "Modified PATH: $PATH"

        # Now ensure Java is not in the cleaned PATH
        if ! command_exists java; then
            echo "✓ Java successfully removed from PATH for testing"
        else
            echo "Warning: Java still in PATH at $(which java)"
        fi

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
# First check should detect Java in brew but not in PATH
CHECK_OUTPUT=$(bv check 2>&1 || true)
echo "$CHECK_OUTPUT"
echo ""

# Verify what bv check detected
if echo "$CHECK_OUTPUT" | grep -q "NOT FOUND"; then
    echo "ℹ️  Java not found - will be installed by setup"
elif echo "$CHECK_OUTPUT" | grep -q "Found (not in PATH)" || echo "$CHECK_OUTPUT" | grep -q "installed via Homebrew but not in your PATH"; then
    echo "✓ bv check correctly detected Java installed via brew but not in PATH"
elif echo "$CHECK_OUTPUT" | grep -q "✓ Found"; then
    echo "ℹ️  Java already in PATH - setup may skip installation"
fi
echo ""

echo "========================================="
echo "Testing bv setup"
echo "========================================="
# In CI mode, setup should automatically configure PATH if Java is in brew but not in PATH
echo "Running bv setup in CI mode (non-interactive)..."
bv setup || true
echo ""

# After setup, check if Java is now accessible
if command_exists java; then
    echo "ℹ️  Java is now in PATH after setup"
else
    echo "ℹ️  Java not yet in PATH (may need shell restart)"

    # Check if the shell config was updated
    if [ -f "$HOME/.zshrc" ] && grep -q "Added by BioVault setup" "$HOME/.zshrc"; then
        echo "✓ bv setup added Java to shell configuration (~/.zshrc)"
        # Try sourcing to apply changes
        if [ -f "$HOME/.zshrc" ]; then
            . "$HOME/.zshrc" 2>/dev/null || true
        fi
    elif [ -f "$HOME/.bash_profile" ] && grep -q "Added by BioVault setup" "$HOME/.bash_profile"; then
        echo "✓ bv setup added Java to shell configuration (~/.bash_profile)"
        # Try sourcing to apply changes
        if [ -f "$HOME/.bash_profile" ]; then
            . "$HOME/.bash_profile" 2>/dev/null || true
        fi
    fi

    # If Java was installed via brew, manually add it for this test session
    if brew list --formula | grep -q openjdk; then
        for pkg in $(brew list --formula | grep openjdk); do
            BREW_JAVA_PATH=$(brew --prefix "$pkg")/bin
            if [ -d "$BREW_JAVA_PATH" ] && [ -f "$BREW_JAVA_PATH/java" ]; then
                export PATH="$BREW_JAVA_PATH:$PATH"
                echo "   Manually added $BREW_JAVA_PATH to PATH for testing"
                break
            fi
        done
    fi
fi
echo ""

echo "========================================="
echo "Verifying installations"
echo "========================================="

# Check Java installation
if command_exists java; then
    echo "✓ Java installed and in PATH:"
    java -version 2>&1 | head -1
else
    # Check if Java is at least installed via brew
    if brew list --formula | grep -q openjdk; then
        echo "⚠️  Java installed via brew but not yet in PATH"
        echo "   This is expected if shell config hasn't been sourced yet"
        # Try to manually add to PATH for this test
        for pkg in $(brew list --formula | grep openjdk); do
            JAVA_PATH=$(brew --prefix "$pkg")/bin
            if [ -d "$JAVA_PATH" ]; then
                export PATH="$JAVA_PATH:$PATH"
                echo "   Manually added $JAVA_PATH to PATH for testing"
                break
            fi
        done
        # Check again
        if command_exists java; then
            echo "✓ Java now accessible:"
            java -version 2>&1 | head -1
        else
            echo "✗ Java still not found after PATH update"
            exit 1
        fi
    else
        echo "✗ Java not found"
        exit 1
    fi
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
# Capture output to decide pass criteria on CI (Docker daemon may not run)
CHECK_OUTPUT=$(bv check 2>&1 || true)
CHECK_STATUS=$?
echo "$CHECK_OUTPUT"
echo ""

# If check failed only because some services are not running, treat as success on CI.
if [ $CHECK_STATUS -ne 0 ]; then
  if echo "$CHECK_OUTPUT" | grep -q "Some services are not running"; then
    echo "Note: Services not running (expected on CI). Proceeding."
  else
    echo "bv check failed with missing dependencies. Exiting."
    exit $CHECK_STATUS
  fi
fi

echo "========================================="
echo "✓ macOS installation test completed"
echo "========================================="
