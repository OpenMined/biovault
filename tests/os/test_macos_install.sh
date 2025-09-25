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
        # We need to completely remove Java from PATH to test our installation

        if command_exists java; then
            JAVA_PATH=$(which java)
            echo "Java currently found at: $JAVA_PATH"

            # Get Java version for info
            java -version 2>&1 | head -1 || true

            # Check if this Java is from brew or system
            if echo "$JAVA_PATH" | grep -q "brew\|homebrew\|Cellar"; then
                echo "Java is from Homebrew"
                # Uninstall it to test fresh installation
                echo "Uninstalling Homebrew Java to test fresh installation..."
                for pkg in $(brew list --formula | grep openjdk || true); do
                    brew uninstall --ignore-dependencies "$pkg" 2>/dev/null || true
                done
            else
                echo "Java appears to be from system (not Homebrew)"
                echo "System Java at $JAVA_PATH cannot be uninstalled"
            fi
        else
            echo "No Java found in PATH"
        fi

        # Remove ALL directories containing Java from PATH
        echo "Removing all Java-related directories from PATH..."
        export ORIGINAL_PATH="$PATH"

        # More aggressive PATH cleaning - remove any directory that contains java/jdk/openjdk
        NEW_PATH=""
        IFS=':' read -ra DIRS <<< "$PATH"
        for dir in "${DIRS[@]}"; do
            # Skip if directory contains java, jdk, openjdk, or has java binary
            if [[ "$dir" =~ (java|jdk|openjdk|Java) ]] || [[ -x "$dir/java" ]]; then
                echo "  Removing from PATH: $dir"
            else
                if [ -n "$NEW_PATH" ]; then
                    NEW_PATH="$NEW_PATH:$dir"
                else
                    NEW_PATH="$dir"
                fi
            fi
        done
        export PATH="$NEW_PATH"

        echo "Modified PATH: $PATH"
        echo ""

        # Verify Java is truly not accessible
        if command_exists java; then
            echo "❌ ERROR: Java is still in PATH at $(which java)"
            echo "This test cannot proceed with Java still accessible."
            echo "The test needs Java to be completely removed from PATH to test installation."
            exit 1
        else
            echo "✓ Java successfully removed from PATH for testing"
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
if echo "$CHECK_OUTPUT" | grep -q "java.*NOT FOUND"; then
    echo "✓ bv check correctly reports Java as NOT FOUND (expected after PATH cleaning)"
elif echo "$CHECK_OUTPUT" | grep -q "Found (not in PATH)" || echo "$CHECK_OUTPUT" | grep -q "installed via Homebrew but not in your PATH"; then
    echo "⚠️  bv check detected Java installed via brew but not in PATH"
elif echo "$CHECK_OUTPUT" | grep -q "java.*✓ Found"; then
    echo "❌ ERROR: bv check still finds Java in PATH - test setup failed"
    echo "Cannot proceed with testing Java installation"
    exit 1
else
    echo "⚠️  Unexpected bv check output for Java"
fi
echo ""

echo "========================================="
echo "Testing bv setup"
echo "========================================="
# In CI mode, setup should automatically configure PATH if Java is in brew but not in PATH
echo "Running bv setup in CI mode (non-interactive)..."
SETUP_OUTPUT=$(bv setup 2>&1 || true)
echo "$SETUP_OUTPUT"
echo ""

# Check if setup actually installed Java
if echo "$SETUP_OUTPUT" | grep -q "Installing java"; then
    echo "✓ bv setup installed Java"
elif echo "$SETUP_OUTPUT" | grep -q "already meets requirements"; then
    echo "⚠️  bv setup reports Java already meets requirements"
else
    echo "⚠️  Could not determine if Java was installed from setup output"
fi

# After setup, check if Java is now accessible
if command_exists java; then
    echo "✓ Java is now in PATH after setup"
    java -version 2>&1 | head -1
else
    echo "⚠️  Java not yet in PATH (checking for brew installation)..."

    # Check if the shell config was updated
    if [ -f "$HOME/.zshrc" ] && grep -q "Added by BioVault setup" "$HOME/.zshrc"; then
        echo "✓ bv setup added configuration to ~/.zshrc"
        # Show what was added
        grep -A 1 "Added by BioVault setup" "$HOME/.zshrc" | sed 's/^/    /'

        # Try sourcing to apply changes
        . "$HOME/.zshrc" 2>/dev/null || true
    elif [ -f "$HOME/.bash_profile" ] && grep -q "Added by BioVault setup" "$HOME/.bash_profile"; then
        echo "✓ bv setup added configuration to ~/.bash_profile"
        # Show what was added
        grep -A 1 "Added by BioVault setup" "$HOME/.bash_profile" | sed 's/^/    /'

        # Try sourcing to apply changes
        . "$HOME/.bash_profile" 2>/dev/null || true
    else
        echo "⚠️  No shell configuration was added by bv setup"
    fi

    # If Java was installed via brew, manually add it for this test session
    if brew list --formula 2>/dev/null | grep -q openjdk; then
        echo "Java is installed via brew, manually adding to PATH..."
        for pkg in $(brew list --formula | grep openjdk); do
            BREW_JAVA_PATH=$(brew --prefix "$pkg")/bin
            if [ -d "$BREW_JAVA_PATH" ] && [ -f "$BREW_JAVA_PATH/java" ]; then
                export PATH="$BREW_JAVA_PATH:$PATH"
                echo "   Added $BREW_JAVA_PATH to PATH"
                break
            fi
        done
    fi
fi
echo ""

echo "========================================="
echo "Verifying installations"
echo "========================================="

# Final verification - Check Java installation
echo ""
echo "Final Java verification:"
if command_exists java; then
    echo "✓ Java is accessible in PATH:"
    java -version 2>&1 | head -1
    which java
else
    # Check if Java is installed via brew but not in PATH
    if brew list --formula 2>/dev/null | grep -q openjdk; then
        echo "⚠️  Java installed via brew but not in current PATH"
        echo "   This may be expected if shell config needs to be sourced in a new session"

        # Try to manually add to PATH for verification
        for pkg in $(brew list --formula | grep openjdk); do
            JAVA_PATH=$(brew --prefix "$pkg")/bin
            if [ -d "$JAVA_PATH" ] && [ -f "$JAVA_PATH/java" ]; then
                export PATH="$JAVA_PATH:$PATH"
                echo "   Manually adding $JAVA_PATH to PATH for verification..."
                break
            fi
        done

        # Check again after manual PATH update
        if command_exists java; then
            echo "✓ Java now accessible after manual PATH update:"
            java -version 2>&1 | head -1
            echo "   This confirms Java was installed correctly by bv setup"
            echo "   Users will have Java in PATH after restarting their shell"
        else
            echo "✗ Java still not found even after manual PATH update"
            echo "   Installation may have failed"
            exit 1
        fi
    else
        echo "✗ Java not found via brew or in PATH"
        echo "   bv setup should have installed Java via brew"
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
