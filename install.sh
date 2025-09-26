#!/bin/bash

# BioVault (bv) installer script
# Usage: curl -sSL https://raw.githubusercontent.com/openmined/biovault/main/install.sh | bash

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Detect OS and architecture
detect_platform() {
    local os
    local arch
    
    # Detect OS
    case "$(uname -s)" in
        Linux*)     os="linux" ;;
        Darwin*)    os="macos" ;;
        CYGWIN*|MINGW*|MSYS*) os="windows" ;;
        *)          
            print_error "Unsupported operating system: $(uname -s)"
            exit 1
            ;;
    esac
    
    # Detect architecture
    case "$(uname -m)" in
        x86_64|amd64)   arch="x86_64" ;;
        aarch64|arm64)  arch="aarch64" ;;
        *)              
            print_error "Unsupported architecture: $(uname -m)"
            exit 1
            ;;
    esac
    
    echo "${os}-${arch}"
}

# Get latest release version
get_latest_version() {
    local api_url="https://api.github.com/repos/openmined/biovault/releases/latest"
    
    # Try to get version from GitHub API
    if command -v curl >/dev/null 2>&1; then
        curl -s "$api_url" | grep '"tag_name":' | sed -E 's/.*"tag_name": "([^"]+)".*/\1/' | head -1
    elif command -v wget >/dev/null 2>&1; then
        wget -qO- "$api_url" | grep '"tag_name":' | sed -E 's/.*"tag_name": "([^"]+)".*/\1/' | head -1
    else
        print_error "Neither curl nor wget is available. Please install one of them."
        exit 1
    fi
}

# Download and install BioVault
install_biovault() {
    local platform="$1"
    local version="$2"
    local install_dir="${3:-/usr/local/bin}"
    
    # Map platform to target architecture
    local target
    case "$platform" in
        linux-x86_64)   target="x86_64-unknown-linux-musl" ;;
        linux-aarch64)  target="aarch64-unknown-linux-musl" ;;
        macos-x86_64)   target="x86_64-apple-darwin" ;;
        macos-aarch64)  target="aarch64-apple-darwin" ;;
        windows-x86_64) target="x86_64-pc-windows-msvc" ;;
        windows-aarch64) target="aarch64-pc-windows-msvc" ;;
        *)
            print_error "Unsupported platform: $platform"
            exit 1
            ;;
    esac
    
    # Construct download URL for the tarball/zip
    local archive_name="bv-${target}"
    local archive_ext="tar.gz"
    if [[ "$platform" == *"windows"* ]]; then
        archive_ext="zip"
    fi
    archive_name="${archive_name}.${archive_ext}"
    echo "archive_name: ${archive_name}"
    
    local download_url="https://github.com/openmined/biovault/releases/download/${version}/${archive_name}"
    echo "download_url: ${download_url}"
    local temp_dir="/tmp/biovault-install-$$"
    echo "temp_dir: ${temp_dir}"
    mkdir -p "$temp_dir"
    
    print_status "Downloading BioVault ${version} for ${platform}..."
    
    # Download the archive
    local temp_archive="${temp_dir}/${archive_name}"
    echo "temp_archive: ${temp_archive}"
    if command -v curl >/dev/null 2>&1; then
        curl -L -o "$temp_archive" "$download_url"
    elif command -v wget >/dev/null 2>&1; then
        wget -O "$temp_archive" "$download_url"
    else
        print_error "Neither curl nor wget is available."
        rm -rf "$temp_dir"
        exit 1
    fi
    
    # Verify download
    if [[ ! -f "$temp_archive" ]]; then
        print_error "Failed to download BioVault archive"
        rm -rf "$temp_dir"
        exit 1
    fi
    
    # Extract the binary
    print_status "Extracting BioVault (bv)..."
    cd "$temp_dir"
    if [[ "$archive_ext" == "tar.gz" ]]; then
        echo "tar -xzf $archive_name"
        if ! tar -xzf "$archive_name"; then
            print_error "Extraction failed: Archive is not in gzip format"
            rm -rf "$temp_dir"
            exit 1
        fi
    elif [[ "$archive_ext" == "zip" ]]; then
        if ! unzip -q "$archive_name"; then
            print_error "Extraction failed: Archive is not in zip format"
            rm -rf "$temp_dir"
            exit 1
        fi
    fi
    
    # Find the binary
    local binary_name="bv"
    if [[ "$platform" == *"windows"* ]]; then
        binary_name="bv.exe"
    fi
    
    if [[ ! -f "$binary_name" ]]; then
        print_error "Binary not found in archive"
        rm -rf "$temp_dir"
        exit 1
    fi
    
    # Make executable
    chmod +x "$binary_name"
    
    # Install to system path
    local target_file="${install_dir}/bv"
    
    print_status "Installing to ${target_file}..."
    
    # Try to install to system directory
    if [[ -w "$install_dir" ]]; then
        mv "$binary_name" "$target_file"
    else
        # Use sudo if directory is not writable
        print_status "Requesting sudo permission to install to ${install_dir}..."
        sudo mv "$binary_name" "$target_file"
    fi
    
    # Cleanup
    cd - >/dev/null
    rm -rf "$temp_dir"
    
    print_success "bv installed successfully!"
}

# Verify installation
verify_installation() {
    local install_dir="$1"

    # Check if binary exists in install dir
    if [[ ! -f "${install_dir}/bv" ]]; then
        print_error "Installation failed. Binary not found at ${install_dir}/bv"
        return 1
    fi

    # Check if bv is directly available in current session
    if command -v bv >/dev/null 2>&1; then
        local installed_version
        installed_version=$(bv --version 2>/dev/null | head -1 || echo "unknown")
        print_success "Installation verified: ${installed_version}"
        print_status ""
        print_status "🎉 BioVault is ready to use!"
        print_status ""
        print_status "Next steps:"
        print_status "  bv check          # Check system dependencies"
        print_status "  bv setup          # Install missing dependencies"
        print_status "  bv init <email>   # Initialize BioVault"
        print_status ""
        print_status "For more information, run: bv --help"
        return 0
    else
        # Binary exists but not in PATH yet
        print_status "bv was installed to ${install_dir}"
        print_status ""
        configure_shell_path "${install_dir}"

        # After configuring, check again if it's available now
        if command -v bv >/dev/null 2>&1; then
            local installed_version
            installed_version=$(bv --version 2>/dev/null | head -1 || echo "unknown")
            print_success "Installation verified: ${installed_version}"
            print_status ""
            print_status "🎉 BioVault is ready to use!"
            print_status ""
            print_status "Next steps:"
            print_status "  bv check          # Check system dependencies"
            print_status "  bv setup          # Install missing dependencies"
            print_status "  bv init <email>   # Initialize BioVault"
        fi
        return 0
    fi
}

# Detect user's shell
detect_shell() {
    local shell_name=""
    if [[ -n "$SHELL" ]]; then
        shell_name=$(basename "$SHELL")
    fi

    # Fallback detection
    if [[ -z "$shell_name" ]]; then
        if [[ -f "/etc/passwd" ]]; then
            shell_name=$(grep "^$(whoami):" /etc/passwd | cut -d: -f7 | xargs basename)
        fi
    fi

    echo "${shell_name:-bash}"
}

# Get shell config file
get_shell_config() {
    local shell_name="$1"
    local config_file=""

    case "$shell_name" in
        zsh)
            config_file="$HOME/.zshrc"
            ;;
        bash)
            # On macOS, use .bash_profile for login shells
            if [[ "$(uname)" == "Darwin" ]]; then
                config_file="$HOME/.bash_profile"
            else
                config_file="$HOME/.bashrc"
            fi
            ;;
        fish)
            config_file="$HOME/.config/fish/config.fish"
            ;;
        *)
            config_file="$HOME/.profile"
            ;;
    esac

    echo "$config_file"
}

# Configure shell PATH
configure_shell_path() {
    local install_dir="$1"
    local shell_name
    local config_file

    shell_name=$(detect_shell)
    config_file=$(get_shell_config "$shell_name")

    print_status "Detected shell: $shell_name"
    print_status "Shell config file: $config_file"

    # Check if we're in CI/non-interactive mode or piped install
    local is_interactive=true
    if [[ -n "$CI" ]] || [[ -n "$GITHUB_ACTIONS" ]] || [[ ! -t 0 ]]; then
        is_interactive=false
    fi

    # For piped installs (curl | bash), we want to auto-configure
    local should_configure=false

    if [[ "$is_interactive" == "false" ]]; then
        # In non-interactive mode, auto-configure for piped installs
        if [[ -z "$CI" ]] && [[ -z "$GITHUB_ACTIONS" ]]; then
            # This is a piped install (curl | bash), not CI
            print_status "Auto-configuring PATH for piped installation..."
            should_configure=true
        else
            # This is CI, just show instructions
            print_status "Running in CI mode."
            print_warning "Please add the following to your shell config manually:"
            if [[ "$shell_name" == "fish" ]]; then
                print_status "  set -gx PATH $install_dir \$PATH"
            else
                print_status "  export PATH=\"$install_dir:\$PATH\""
            fi
            return
        fi
    else
        # Interactive mode - ask user
        print_status ""
        print_status "Would you like to automatically add $install_dir to your PATH? [Y/n]: "
        read -r response
        response=${response:-Y}
        if [[ "$response" =~ ^[Yy] ]]; then
            should_configure=true
        fi
    fi

    if [[ "$should_configure" == "true" ]]; then
        # Create config file if it doesn't exist
        if [[ ! -f "$config_file" ]]; then
            touch "$config_file"
        fi

        # Check if PATH update already exists
        if grep -q "$install_dir" "$config_file" 2>/dev/null; then
            print_status "PATH configuration already exists in $config_file"
        else
            # Add PATH configuration
            echo "" >> "$config_file"
            echo "# Added by BioVault installer on $(date)" >> "$config_file"

            if [[ "$shell_name" == "fish" ]]; then
                echo "set -gx PATH $install_dir \$PATH" >> "$config_file"
            else
                echo "export PATH=\"$install_dir:\$PATH\"" >> "$config_file"
            fi

            print_success "PATH configuration added to $config_file"
        fi

        # Export PATH for current session
        export PATH="$install_dir:$PATH"
        print_success "PATH updated for current session"

        print_status ""
        print_success "Installation complete! You can now run 'bv' commands."
        print_status "Note: New terminal sessions will have 'bv' in PATH automatically."
    else
        print_status "Skipping PATH configuration."
        print_status "To add manually, add the following to $config_file:"
        if [[ "$shell_name" == "fish" ]]; then
            print_status "  set -gx PATH $install_dir \$PATH"
        else
            print_status "  export PATH=\"$install_dir:\$PATH\""
        fi
        print_status "Then run: source $config_file"
    fi
}

# Check prerequisites
check_prerequisites() {
    print_status "Checking prerequisites..."

    # Currently no specific prerequisites required
    print_success "All prerequisites met!"
}

# Main installation function
main() {
    print_status "BioVault (bv) installer"
    print_status "======================="
    
    # Parse command line arguments
    local custom_install_dir=""
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --prefix|--install-dir)
                custom_install_dir="$2"
                shift 2
                ;;
            --help|-h)
                echo "Usage: $0 [OPTIONS]"
                echo "Options:"
                echo "  --prefix, --install-dir <DIR>  Install to specific directory"
                echo "  --help, -h                     Show this help message"
                echo ""
                echo "Examples:"
                echo "  $0                             # Install to default location"
                echo "  $0 --prefix ~/.local/bin       # Install to ~/.local/bin"
                echo "  $0 --install-dir ~/bin         # Install to ~/bin"
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done
    
    # check_prerequisites
    
    # Detect platform
    local platform
    platform=$(detect_platform)
    print_status "Detected platform: ${platform}"
    
    # Get latest version
    local version
    version=$(get_latest_version)
    if [[ -z "$version" ]]; then
        print_error "Failed to get latest version information"
        exit 1
    fi
    print_status "Latest version: ${version}"
    
    # Determine install directory - prefer user directories to avoid sudo
    local install_dir=""
    
    # Use custom directory if specified
    if [[ -n "$custom_install_dir" ]]; then
        install_dir="$custom_install_dir"
        if [[ ! -d "$install_dir" ]]; then
            print_status "Creating directory: $install_dir"
            mkdir -p "$install_dir"
        fi
        print_status "Using specified directory: $install_dir"
    else
        # First, try to find a writable directory that's already in PATH
        local path_dirs
        IFS=':' read -ra path_dirs <<< "$PATH"

        # Priority 1: Check directories already in PATH that are writable
        # Prefer user directories to avoid needing sudo
        local checked_dirs=()

        for dir in "${path_dirs[@]}"; do
            # Skip empty entries and duplicates
            [[ -z "$dir" ]] && continue
            [[ " ${checked_dirs[@]} " =~ " ${dir} " ]] && continue
            checked_dirs+=("$dir")

            # Check if directory exists and is writable
            if [[ -d "$dir" ]] && [[ -w "$dir" ]]; then
                # Prefer user directories over system directories
                if [[ "$dir" == "$HOME"* ]]; then
                    install_dir="$dir"
                    print_status "Found writable user directory in PATH: $install_dir"
                    break
                elif [[ -z "$install_dir" ]]; then
                    # Keep system directory as fallback
                    install_dir="$dir"
                fi
            fi
        done

        # If we found a system directory but not a user directory, use it
        if [[ -n "$install_dir" ]] && [[ "$install_dir" != "$HOME"* ]]; then
            print_status "Found writable system directory in PATH: $install_dir"
        fi

        # Priority 2: Check common general-purpose bin directories that might be in PATH
        # Only use general bin directories, not package-manager-specific ones
        if [[ -z "$install_dir" ]]; then
            local common_user_dirs=(
                "$HOME/.local/bin"      # XDG standard (Linux/macOS)
                "$HOME/bin"             # Traditional Unix user bin
                "$HOME/.bin"            # Hidden user bin (some distros)
            )

            for dir in "${common_user_dirs[@]}"; do
                # Check if it's already in PATH
                if [[ ":$PATH:" == *":$dir:"* ]]; then
                    # Create if needed (user directory)
                    if [[ ! -d "$dir" ]]; then
                        print_status "Creating directory: $dir"
                        mkdir -p "$dir"
                    fi
                    # Check if writable
                    if [[ -w "$dir" ]]; then
                        install_dir="$dir"
                        print_status "Using existing PATH directory: $install_dir"
                        break
                    fi
                fi
            done
        fi

        # Priority 3: Try system directories if no user directory found
        if [[ -z "$install_dir" ]]; then
            local system_dirs=("/usr/local/bin" "/opt/bin")
            for dir in "${system_dirs[@]}"; do
                if [[ ":$PATH:" == *":$dir:"* ]] && [[ -d "$dir" ]] && [[ -w "$dir" ]]; then
                    install_dir="$dir"
                    print_status "Using system directory: $install_dir (already in PATH)"
                    break
                fi
            done
        fi
        
        # If no directory in PATH is writable, try to use ~/.local/bin
        if [[ -z "$install_dir" ]]; then
            install_dir="$HOME/.local/bin"
            if [[ ! -d "$install_dir" ]]; then
                print_status "Creating local bin directory: $install_dir"
                mkdir -p "$install_dir"
            fi

            # Note: We'll handle PATH configuration after installation
            if [[ ":$PATH:" != *":$install_dir:"* ]]; then
                print_warning "Directory $install_dir is not currently in PATH"
                print_status "We'll help you configure this after installation."
            fi
        fi
        
        # Fall back to system directory with sudo if needed
        if [[ ! -w "$install_dir" ]]; then
            install_dir="/usr/local/bin"
            print_warning "Will need sudo to install to $install_dir"
        fi
    fi
    
    # Install biovault
    install_biovault "$platform" "$version" "$install_dir"

    # Verify installation
    verify_installation "$install_dir"
}

# Run main function
main "$@"