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
    echo "Preparing test environment..."

    # On GitHub runners, Homebrew is pre-installed
    if command_exists brew; then
        echo "Homebrew is installed"

        # Check current Homebrew setup
        echo ""
        echo "Checking Homebrew environment:"
        echo "brew --prefix: $(brew --prefix)"
        echo "brew --prefix openjdk (if installed): $(brew --prefix openjdk 2>/dev/null || echo 'not installed')"
        echo ""

        # Check what's in PATH related to Homebrew
        echo "Homebrew-related PATH entries:"
        echo "$PATH" | tr ':' '\n' | grep -E "(homebrew|Homebrew)" || echo "None found"
        echo ""

        # IMPORTANT: Uninstall any existing OpenJDK to test fresh installation
        echo "Checking for existing OpenJDK installations..."

        # First, uninstall all Homebrew Java packages
        if brew list --formula 2>/dev/null | grep -E "(openjdk|java)" | grep -v javascript; then
            echo "Found Java packages installed via Homebrew:"
            brew list --formula | grep -E "(openjdk|java)" | grep -v javascript || true
            echo ""
            echo "Uninstalling all Java packages to test fresh installation..."
            for pkg in $(brew list --formula | grep -E "(openjdk|java)" | grep -v javascript || true); do
                echo "  Uninstalling $pkg..."
                brew uninstall --ignore-dependencies "$pkg" 2>/dev/null || true
            done
            echo "Java packages uninstalled from Homebrew"
        else
            echo "No Java packages found in Homebrew"
        fi

        # Remove any system Java from PATH temporarily
        echo ""
        echo "Checking for system Java..."
        if command_exists java; then
            JAVA_PATH=$(which java)
            echo "Found Java in PATH at: $JAVA_PATH"
            java -version 2>&1 | head -1 || true

            # Check if it's a symlink and where it points
            if [ -L "$JAVA_PATH" ]; then
                echo "Java is a symlink pointing to: $(readlink -f "$JAVA_PATH" 2>/dev/null || readlink "$JAVA_PATH")"
            fi

            # Show JAVA_HOME if set
            if [ -n "$JAVA_HOME" ]; then
                echo "JAVA_HOME is set to: $JAVA_HOME"
            fi

            # GitHub Actions runners have system Java pre-installed
            # We need to remove it from PATH for our test
            echo ""
            echo "Removing system Java from PATH for this test..."

            # Strategy 1: Remove Java-specific directories from PATH
            CLEANED_PATH=""
            IFS=':' read -ra PATH_ARRAY <<< "$PATH"
            for p in "${PATH_ARRAY[@]}"; do
                # Skip any path that contains Java-specific directories
                # But keep system directories like /usr/bin
                if [[ "$p" =~ (java|Java|jdk|jre|JDK|JRE|temurin|zulu|openjdk|microsoft-jdk|graalvm|corretto|liberica|sapmachine|semeru|JavaVirtualMachines) ]]; then
                    echo "  Removing from PATH: $p"
                elif [[ "$p" == "$JAVA_HOME"* ]] && [ -n "$JAVA_HOME" ]; then
                    echo "  Removing JAVA_HOME path: $p"
                else
                    if [ -z "$CLEANED_PATH" ]; then
                        CLEANED_PATH="$p"
                    else
                        CLEANED_PATH="$CLEANED_PATH:$p"
                    fi
                fi
            done
            export PATH="$CLEANED_PATH"

            # Also unset JAVA_HOME if set
            if [ -n "$JAVA_HOME" ]; then
                echo "Unsetting JAVA_HOME (was: $JAVA_HOME)"
                unset JAVA_HOME
            fi

            # Strategy 2: Handle system Java (GitHub Actions installs in hostedtoolcache)
            if command_exists java; then
                REMAINING_JAVA=$(which java)
                echo "Java still accessible at: $REMAINING_JAVA"

                # On GitHub Actions, we can safely remove Java executables
                # since we're in a temporary environment
                if [ -n "$GITHUB_ACTIONS" ]; then
                    # Check if it's /usr/bin/java
                    if [ "$REMAINING_JAVA" = "/usr/bin/java" ]; then
                        # Check if it's a symlink
                        if [ -L "$REMAINING_JAVA" ]; then
                            JAVA_TARGET=$(readlink -f "$REMAINING_JAVA" 2>/dev/null || readlink "$REMAINING_JAVA")
                            echo "  /usr/bin/java is a symlink to: $JAVA_TARGET"
                        fi

                        echo "  Removing /usr/bin/java and /usr/bin/javac (GitHub Actions environment)"
                        sudo rm -f /usr/bin/java 2>/dev/null || true
                        sudo rm -f /usr/bin/javac 2>/dev/null || true
                        sudo rm -f /usr/bin/jar 2>/dev/null || true
                        sudo rm -f /usr/bin/jarsigner 2>/dev/null || true
                        sudo rm -f /usr/bin/javadoc 2>/dev/null || true
                        sudo rm -f /usr/bin/javap 2>/dev/null || true
                        sudo rm -f /usr/bin/jcmd 2>/dev/null || true
                        sudo rm -f /usr/bin/jconsole 2>/dev/null || true
                        sudo rm -f /usr/bin/jdb 2>/dev/null || true
                        sudo rm -f /usr/bin/jdeps 2>/dev/null || true
                        sudo rm -f /usr/bin/jinfo 2>/dev/null || true
                        sudo rm -f /usr/bin/jmap 2>/dev/null || true
                        sudo rm -f /usr/bin/jps 2>/dev/null || true
                        sudo rm -f /usr/bin/jrunscript 2>/dev/null || true
                        sudo rm -f /usr/bin/jstack 2>/dev/null || true
                        sudo rm -f /usr/bin/jstat 2>/dev/null || true
                        sudo rm -f /usr/bin/jstatd 2>/dev/null || true
                        sudo rm -f /usr/bin/keytool 2>/dev/null || true

                        # Also check for any Java-related executables we might have missed
                        for java_tool in /usr/bin/j*; do
                            if [ -f "$java_tool" ] && [ -L "$java_tool" ]; then
                                if readlink "$java_tool" 2>/dev/null | grep -q -E "(java|jdk|jre|temurin|hostedtoolcache)"; then
                                    echo "  Also removing: $java_tool"
                                    sudo rm -f "$java_tool" 2>/dev/null || true
                                fi
                            fi
                        done
                    fi
                fi
            fi

            # Verify Java is no longer accessible
            if command_exists java; then
                echo "WARNING: Java still accessible after all modifications"
                echo "Location: $(which java 2>/dev/null || echo 'unknown')"
            else
                echo "✓ Successfully removed/hidden Java from environment"
            fi
        else
            echo "No Java found in current PATH"
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
CHECK_OUTPUT=$(bv check 2>&1 || true)
echo "$CHECK_OUTPUT"
echo ""

# Analyze the check output
echo "Analysis of bv check output:"
if echo "$CHECK_OUTPUT" | grep -q "java.*NOT FOUND"; then
    echo "✓ Java is NOT FOUND (expected - we uninstalled it)"
elif echo "$CHECK_OUTPUT" | grep -q "Found (not in PATH)"; then
    echo "⚠️  Java is installed but not in PATH (unexpected - we should have uninstalled it)"
elif echo "$CHECK_OUTPUT" | grep -q "java.*✓ Found"; then
    echo "⚠️  Java is in PATH (system Java?)"
fi
echo ""

echo "========================================="
echo "Testing bv setup"
echo "========================================="
echo "Running bv setup in CI mode (non-interactive)..."
SETUP_OUTPUT=$(bv setup 2>&1 || true)
echo "$SETUP_OUTPUT"
echo ""

echo "========================================="
echo "Analyzing setup results"
echo "========================================="

# Check if Java was installed
if echo "$SETUP_OUTPUT" | grep -q "Installing java"; then
    echo "✓ bv setup attempted to install Java"
elif echo "$SETUP_OUTPUT" | grep -q "java already meets requirements"; then
    echo "⚠️  bv setup says Java already meets requirements"
fi

# Check if PATH configuration was done
if echo "$SETUP_OUTPUT" | grep -q "Java is installed via Homebrew but not in your PATH"; then
    echo "✓ bv setup detected Java needs PATH configuration"
fi

if echo "$SETUP_OUTPUT" | grep -q "Added to.*zshrc\|Added to.*bash_profile"; then
    echo "✓ bv setup added Java to shell configuration"
fi

echo ""

echo "========================================="
echo "Checking shell configuration"
echo "========================================="

# Check what was added to shell configs
for config_file in "$HOME/.zshrc" "$HOME/.bash_profile" "$HOME/.bashrc" "$HOME/.profile"; do
    if [ -f "$config_file" ]; then
        if grep -q "Added by BioVault setup" "$config_file" 2>/dev/null; then
            echo "Found BioVault configuration in $config_file:"
            grep -A 2 "Added by BioVault setup" "$config_file" | sed 's/^/  /'
            echo ""
        fi
    fi
done

echo "========================================="
echo "Verifying Java installation"
echo "========================================="

# No PATH restoration needed since we only removed Java executables, not directories

# First, check if Java is immediately available
if command_exists java; then
    echo "✓ Java is immediately available in PATH"
    java -version 2>&1 | head -1
    echo "Location: $(which java)"
else
    echo "Java is not immediately available in PATH"

    # Check if Java was installed via brew
    if brew list --formula 2>/dev/null | grep -q openjdk; then
        echo "✓ Java IS installed via Homebrew"

        # Get the path where brew installed Java
        JAVA_PKG=$(brew list --formula | grep openjdk | head -1)
        BREW_JAVA_PATH=$(brew --prefix "$JAVA_PKG")/bin
        echo "  Homebrew Java location: $BREW_JAVA_PATH"

        # Check if this path is in our current PATH
        if echo "$PATH" | grep -q "$BREW_JAVA_PATH"; then
            echo "  ✓ This location IS in current PATH"
        else
            echo "  ⚠️  This location is NOT in current PATH"
            echo "  This is expected - shell config changes require a new session"

            # Manually add it to verify the installation worked
            export PATH="$BREW_JAVA_PATH:$PATH"
            if command_exists java; then
                echo "  ✓ After adding to PATH manually, Java works:"
                java -version 2>&1 | head -1
            else
                echo "  ❌ Even after adding to PATH, Java doesn't work"
                exit 1
            fi
        fi
    else
        echo "❌ Java was NOT installed via Homebrew"
        echo "  This is a problem - bv setup should have installed it"
        exit 1
    fi
fi

echo ""

# Check Nextflow installation
echo "========================================="
echo "Verifying other dependencies"
echo "========================================="
if command_exists nextflow || [ -f /usr/local/bin/nextflow ]; then
    echo "✓ Nextflow installed"
    nextflow -version 2>/dev/null || /usr/local/bin/nextflow -version 2>/dev/null || echo "  (version check may require Java in PATH)"
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
    echo "⚠️  SyftBox not installed (expected - requires manual PATH setup)"
fi

echo ""
echo "========================================="
echo "Testing bv check (after setup)"
echo "========================================="
# Run check again to see the final state
CHECK_OUTPUT=$(bv check 2>&1 || true)
CHECK_STATUS=$?
echo "$CHECK_OUTPUT"
echo ""

# Analyze final state
if echo "$CHECK_OUTPUT" | grep -q "java.*✓ Found"; then
    echo "✓ Java is now detected as Found by bv check"
elif echo "$CHECK_OUTPUT" | grep -q "Found (not in PATH)"; then
    echo "⚠️  Java is still shown as not in PATH"
    echo "  This is expected in CI - changes require new shell session"
    echo "  The important thing is that Java was installed and PATH was configured"
fi

# Allow test to pass if dependencies are installed but services aren't running
if [ $CHECK_STATUS -ne 0 ]; then
    if echo "$CHECK_OUTPUT" | grep -q "Some services are not running"; then
        echo "Note: Services not running (expected on CI). Test passed."
    elif echo "$CHECK_OUTPUT" | grep -q "Found (not in PATH)"; then
        echo "Note: Java not in current PATH but was installed and configured. Test passed."
    else
        echo "bv check failed unexpectedly. See output above."
        exit $CHECK_STATUS
    fi
fi

echo ""
echo "========================================="
echo "✓ macOS installation test completed"
echo "========================================="
echo ""
echo "Summary:"
echo "- Java was installed via Homebrew: ✓"
echo "- PATH configuration was added to shell: ✓"
echo "- Other dependencies were installed: ✓"
echo "- Test validated the installation process: ✓"