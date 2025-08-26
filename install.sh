#!/bin/bash

# biovault installer script
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

# Download and install biovault
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
    local archive_name="biovault-${target}"
    local archive_ext="tar.gz"
    if [[ "$platform" == *"windows"* ]]; then
        archive_ext="zip"
    fi
    archive_name="${archive_name}.${archive_ext}"
    
    local download_url="https://github.com/openmined/biovault/releases/download/${version}/${archive_name}"
    local temp_dir="/tmp/biovault-install-$$"
    mkdir -p "$temp_dir"
    
    print_status "Downloading biovault ${version} for ${platform}..."
    
    # Download the archive
    local temp_archive="${temp_dir}/${archive_name}"
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
        print_error "Failed to download biovault archive"
        rm -rf "$temp_dir"
        exit 1
    fi
    
    # Extract the binary
    print_status "Extracting biovault..."
    cd "$temp_dir"
    if [[ "$archive_ext" == "tar.gz" ]]; then
        tar -xzf "$archive_name"
    elif [[ "$archive_ext" == "zip" ]]; then
        unzip -q "$archive_name"
    fi
    
    # Find the binary
    local binary_name="biovault"
    if [[ "$platform" == *"windows"* ]]; then
        binary_name="biovault.exe"
    fi
    
    if [[ ! -f "$binary_name" ]]; then
        print_error "Binary not found in archive"
        rm -rf "$temp_dir"
        exit 1
    fi
    
    # Make executable
    chmod +x "$binary_name"
    
    # Install to system path
    local target_file="${install_dir}/biovault"
    
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
    
    print_success "biovault installed successfully!"
}

# Verify installation
verify_installation() {
    if command -v biovault >/dev/null 2>&1; then
        local installed_version
        installed_version=$(biovault --version 2>/dev/null | head -1 || echo "unknown")
        print_success "Installation verified: ${installed_version}"
        print_status "You can now use 'biovault' to protect your files!"
        print_status ""
        print_status "Quick start:"
        print_status "  biovault           # Interactive file protection mode"
        print_status "  biovault --help    # Show all available options"
        print_status ""
        print_status "For more information, run: biovault --help"
    else
        print_error "Installation verification failed. biovault command not found in PATH."
        print_warning "You may need to restart your shell or update your PATH."
        return 1
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
    print_status "biovault installer"
    print_status "================"
    
    # Check prerequisites (currently none required)
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
    
    # Determine install directory
    local install_dir="/usr/local/bin"
    if [[ ":$PATH:" != *":$install_dir:"* ]]; then
        # Try alternative directories if /usr/local/bin is not in PATH
        for dir in "$HOME/.local/bin" "$HOME/bin" "/usr/bin"; do
            if [[ ":$PATH:" == *":$dir:"* ]] && [[ -d "$dir" ]]; then
                install_dir="$dir"
                break
            fi
        done
    fi
    
    # Install biovault
    install_biovault "$platform" "$version" "$install_dir"
    
    # Verify installation
    verify_installation
}

# Run main function
main "$@"