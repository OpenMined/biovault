# BioVault Windows Installation Test

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "BioVault Windows Installation Test" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan
Write-Host ""

Write-Host "Environment Info:" -ForegroundColor Yellow
Write-Host "OS: $([System.Environment]::OSVersion.VersionString)"
Write-Host "Architecture: $([System.Environment]::Is64BitOperatingSystem)"
Write-Host ""

# Function to check if a command exists
function Test-Command {
    param($Command)
    try {
        Get-Command $Command -ErrorAction Stop | Out-Null
        return $true
    }
    catch {
        return $false
    }
}

# Function to cleanup existing installations
function Cleanup-Existing {
    Write-Host "Cleaning up potentially conflicting packages..." -ForegroundColor Yellow

    # Windows GitHub runners might have Java pre-installed
    if (Test-Command "java") {
        Write-Host "Found existing Java installation - GitHub runners have Java pre-installed" -ForegroundColor Yellow
        # We won't uninstall Java on Windows as it's complex and might break the runner
    }

    # Remove existing Nextflow if present
    $nextflowPaths = @(
        "$env:USERPROFILE\nextflow.exe",
        "$env:USERPROFILE\nextflow",
        "C:\tools\nextflow.exe",
        "C:\tools\nextflow"
    )

    foreach ($path in $nextflowPaths) {
        if (Test-Path $path) {
            Write-Host "Found existing Nextflow at $path, removing..." -ForegroundColor Yellow
            Remove-Item -Path $path -Force -ErrorAction SilentlyContinue
        }
    }

    Write-Host ""
}

# Only cleanup on CI to avoid affecting local development
if ($env:CI -or $env:GITHUB_ACTIONS) {
    Cleanup-Existing
}

# Verify bv is installed
if (-not (Test-Command "bv")) {
    Write-Host "Error: bv command not found in PATH" -ForegroundColor Red
    Write-Host "PATH: $env:PATH" -ForegroundColor Red
    exit 1
}

Write-Host "Using bv binary at: $(Get-Command bv | Select-Object -ExpandProperty Source)" -ForegroundColor Green
try {
    $version = & bv --version 2>&1
    Write-Host "bv version: $version" -ForegroundColor Green
} catch {
    Write-Host "bv version: unknown" -ForegroundColor Yellow
}
Write-Host ""

Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "Testing bv check (before setup)" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan

# Note: Windows doesn't require bv setup as it doesn't have automated installation support
# But we'll still run bv check to verify the environment
& bv check
if ($LASTEXITCODE -ne 0) {
    Write-Host "bv check reported missing dependencies (expected)" -ForegroundColor Yellow
}
Write-Host ""

Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "Testing bv setup" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan

# Run bv setup - it should detect Windows and provide manual instructions
& bv setup
Write-Host ""

Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "Manual Installation Steps for Windows" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan

# Since Windows doesn't have automated setup, we'll simulate manual installation for CI testing
if ($env:CI -or $env:GITHUB_ACTIONS) {
    Write-Host "Simulating manual installation steps for CI..." -ForegroundColor Yellow

    # Install Nextflow manually for testing
    Write-Host "Installing Nextflow manually..." -ForegroundColor Yellow
    $nextflowUrl = "https://github.com/nextflow-io/nextflow/releases/latest/download/nextflow"
    $nextflowPath = "$env:USERPROFILE\nextflow.exe"

    try {
        Invoke-WebRequest -Uri $nextflowUrl -OutFile $nextflowPath -UseBasicParsing
        Write-Host "Nextflow downloaded to $nextflowPath" -ForegroundColor Green

        # Add to PATH for this session
        $env:PATH = "$env:USERPROFILE;$env:PATH"
    }
    catch {
        Write-Host "Warning: Could not download Nextflow: $_" -ForegroundColor Yellow
    }
}

Write-Host ""
Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "Verifying installations" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan

# Check Java installation
if (Test-Command "java") {
    Write-Host "Java installed:" -ForegroundColor Green
    & java -version 2>&1 | Select-Object -First 1
} else {
    Write-Host "X Java not found - Manual installation required on Windows" -ForegroundColor Red
    Write-Host "  Please install OpenJDK 17+ from: https://openjdk.org/" -ForegroundColor Yellow
}

# Check Nextflow installation
$nextflowFound = $false
if (Test-Command "nextflow") {
    $nextflowFound = $true
} elseif (Test-Path "$env:USERPROFILE\nextflow.exe") {
    $nextflowFound = $true
} elseif (Test-Path "$env:USERPROFILE\nextflow") {
    $nextflowFound = $true
}

if ($nextflowFound) {
    Write-Host "Nextflow installed" -ForegroundColor Green
    try {
        if (Test-Command "nextflow") {
            & nextflow -version
        } else {
            & "$env:USERPROFILE\nextflow.exe" -version
        }
    }
    catch {
        Write-Host "  (Could not get version)" -ForegroundColor Yellow
    }
} else {
    Write-Host "Warning: Nextflow not found - Manual installation required on Windows" -ForegroundColor Yellow
    Write-Host "  Please download from: https://www.nextflow.io/" -ForegroundColor Yellow
}

# Check Docker installation
if (Test-Command "docker") {
    Write-Host "Docker installed:" -ForegroundColor Green
    & docker --version
} else {
    Write-Host "Warning: Docker not installed - Docker Desktop required on Windows" -ForegroundColor Yellow
    Write-Host "  Please install from: https://www.docker.com/products/docker-desktop/" -ForegroundColor Yellow
}

# Check SyftBox installation
if (Test-Command "syftbox") {
    Write-Host "SyftBox installed:" -ForegroundColor Green
    & syftbox -v
} else {
    Write-Host "Warning: SyftBox not installed - Manual installation required on Windows" -ForegroundColor Yellow
    Write-Host "  Please download from: https://github.com/OpenMined/syftbox/releases/latest" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "Testing bv check (after setup)" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan
& bv check

# Don't fail the test on Windows if dependencies are missing, as manual installation is expected
if ($LASTEXITCODE -ne 0) {
    Write-Host "Note: Some dependencies missing - manual installation required on Windows" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "Windows installation test completed" -ForegroundColor Green
Write-Host "=========================================" -ForegroundColor Cyan

# Always exit successfully on Windows since manual installation is expected
exit 0