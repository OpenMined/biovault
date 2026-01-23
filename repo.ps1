#!/usr/bin/env pwsh
# repo.ps1 - Workspace dependency management tool for Windows
# Version: 1.0.0
# PowerShell equivalent of ./repo bash script
# Keep in sync across: biovault-desktop, syftbox-sdk, biovault, biovault-beaver

# Parse arguments manually to match bash script style (--init, --https, etc.)
$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ManifestFile = if ($env:MANIFEST_FILE) { $env:MANIFEST_FILE } else { "manifest.xml" }

function Resolve-WorkspaceRoot {
    param(
        [string]$StartDir,
        [string]$ManifestName
    )

    $current = $StartDir
    while ($true) {
        if (Test-Path (Join-Path $current $ManifestName)) {
            return $current
        }
        $parent = Split-Path -Parent $current
        if (-not $parent -or $parent -eq $current) {
            return $StartDir
        }
        $current = $parent
    }
}

# Allow running from nested folders (e.g., biovault/) by locating manifest.xml.
$RootDir = Resolve-WorkspaceRoot -StartDir $ScriptDir -ManifestName $ManifestFile
$ManifestPath = Join-Path $RootDir $ManifestFile

function Show-Usage {
    Write-Host @"
Usage: .\repo.ps1 [OPTIONS]

Options:
  (none)              Show repo tree with branch/dirty status
  --init [--https]    Initialize workspace and sync deps (clones repos)
  sync                Re-sync workspace to manifest
  fetch               Fetch remotes for all repos and show ahead/behind
  pull [--rebase]     Pull updates for all repos on branches
  lint [--force] [--verbose]
                      Run lint.sh in dirty repos (parallel, quiet on success)
                      --force runs lint.sh in all repos that have it
                      --verbose prints lint output even on success
  test [--force]
                      Run fast unit tests in dirty repos (quiet on success)
                      --force runs tests in all repos with fast tests
  ssh                 Rewrite remotes to SSH for all repos
  main                Checkout main in all repos (no reset)
  switch [-b] <branch> <targets...>
                      Checkout a branch across selected repos
  track <target> <rev>
                      Update manifest revision for a repo (branch/tag/sha)
  checkout [--reset] <rev> <target>
                      Checkout a revision in a repo or all repos
                      target: repo path, repo name, or "all"
  tools [--install] [--firewall]
                      Check required tools (--install to install missing)
                      --firewall adds Windows Firewall allow rules for Node/MinIO
  --help, -h          Show this help

Examples:
  .\repo.ps1                    # Show workspace status
  .\repo.ps1 --init             # Initialize workspace (clone all repos)
  .\repo.ps1 --init --https     # Initialize using HTTPS URLs
  .\repo.ps1 tools              # Check development tools
  .\repo.ps1 tools --firewall   # Add Windows Firewall rules for Node/MinIO
  .\repo.ps1 fetch              # Fetch all repos
  .\repo.ps1 pull               # Pull all repos
  .\repo.ps1 switch -b feature all  # Create/switch to branch in all repos
"@
}

# Tool definitions: name, check_command, winget_id, description, required
# Based on CI requirements in .github/workflows/test.yml
$RequiredTools = @(
    # Core build tools (required)
    @{ Name = "git"; Check = "git --version"; Winget = "Git.Git"; Desc = "Git version control"; Required = $true }
    # MSVC Build Tools - required for Rust to compile native code on Windows
    # Note: cl.exe isn't in PATH normally, but Rust finds it via VS integration
    @{ Name = "msvc"; Check = "if (Test-Path 'C:\Program Files*\Microsoft Visual Studio\2022\BuildTools\VC\Tools\MSVC') { 'VS 2022 Build Tools' } else { throw 'not found' }"; Winget = "Microsoft.VisualStudio.2022.BuildTools"; WingetArgs = "--override `"--add Microsoft.VisualStudio.Workload.VCTools --includeRecommended --passive`""; Desc = "MSVC Build Tools (for Rust)"; Required = $true }
    @{ Name = "rustc"; Check = "rustc --version"; Winget = "Rustlang.Rust.MSVC"; Desc = "Rust compiler for Tauri"; Required = $true }
    @{ Name = "cargo"; Check = "cargo --version"; Winget = "Rustlang.Rust.MSVC"; Desc = "Rust package manager"; Required = $true }
    @{ Name = "go"; Check = "go version"; Winget = $null; ManualInstall = "https://go.dev/dl/go1.24.3.windows-amd64.msi"; RequiredVersion = "1.24"; Desc = "Go 1.24.x for syftbox (1.25+ has issues)"; Required = $true }
    @{ Name = "node"; Check = "node --version"; Winget = "OpenJS.NodeJS.LTS"; Desc = "Node.js for frontend"; Required = $true }
    @{ Name = "npm"; Check = "npm --version"; Winget = "OpenJS.NodeJS.LTS"; Desc = "Node package manager"; Required = $true }
    @{ Name = "webview2"; Check = "if ((Get-ChildItem -Path 'C:\\Program Files (x86)\\Microsoft\\EdgeWebView\\Application' -Filter 'msedgewebview2.exe' -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1) -or (Get-ChildItem -Path 'C:\\Program Files\\Microsoft\\EdgeWebView\\Application' -Filter 'msedgewebview2.exe' -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1)) { 'WebView2 Runtime' } else { throw 'not found' }"; Winget = "Microsoft.EdgeWebView2Runtime"; Desc = "WebView2 runtime for Tauri"; Required = $true }
    @{ Name = "python"; Check = "python --version"; Winget = $null; UseUv = $true; Desc = "Python 3.11 for scripts"; Required = $true }
    @{ Name = "pip"; Check = "python -m pip --version"; Winget = $null; UseUv = $true; Desc = "Python package manager"; Required = $true }
    @{ Name = "protoc"; Check = "protoc --version"; Winget = "Google.Protobuf"; Desc = "Protocol Buffers compiler"; Required = $true }
    # Testing tools (required for E2E tests)
    @{ Name = "bun"; Check = "bun --version"; Winget = "Oven-sh.Bun"; Desc = "Fast JS runtime for Playwright"; Required = $true }
    # Pipeline tools (required for pipelines-collab)
    @{ Name = "docker"; Check = "docker --version"; Winget = "Docker.DockerDesktop"; Desc = "Docker for Nextflow pipelines"; Required = $true }
    @{ Name = "psexec"; Check = "if (Get-Command psexec -ErrorAction SilentlyContinue) { 'PsExec' } else { throw 'not found' }"; Winget = $null; ManualInstall = "https://learn.microsoft.com/sysinternals/downloads/psexec"; Desc = "Sysinternals PsExec for launching GUI apps from SSH"; Required = $false }
    # Optional but recommended
    @{ Name = "uv"; Check = "uv --version"; Winget = "astral-sh.uv"; Desc = "Fast Python package manager"; Required = $false }
)

# Cargo-installed tools (installed via cargo install, not winget)
$CargoTools = @(
    @{ Name = "bvs"; Check = "bvs --version"; Crate = "biosynth"; Desc = "Biosynth CLI for synthetic data" }
)

# Ensure tool directories are in PATH for this session
$toolDirs = @(
    "$env:USERPROFILE\.cargo\bin"
    "$env:USERPROFILE\.local\bin"
    "$env:PROGRAMFILES\Go\bin"
    "$env:USERPROFILE\go\bin"
    "$env:LOCALAPPDATA\Microsoft\WinGet\Links"
    "$env:USERPROFILE\.bun\bin"
)
foreach ($dir in $toolDirs) {
    if ((Test-Path $dir) -and ($env:PATH -notlike "*$dir*")) {
        $env:PATH = "$dir;$env:PATH"
    }
}

# Add uv-managed Python to PATH if available
$uvPythonBase = "$env:APPDATA\uv\python"
if (Test-Path $uvPythonBase) {
    $pythonDir = Get-ChildItem -Path $uvPythonBase -Directory -Filter "cpython-3.11*" | Select-Object -First 1
    if ($pythonDir -and ($env:PATH -notlike "*$($pythonDir.FullName)*")) {
        $env:PATH = "$($pythonDir.FullName);$($pythonDir.FullName)\Scripts;$env:PATH"
    }
}

function Remove-WindowsStorePythonStubs {
    # Windows Store installs stub executables that intercept python/python3 commands
    # and redirect to the Microsoft Store. These need to be removed for real Python to work.
    $stubs = @(
        "$env:LOCALAPPDATA\Microsoft\WindowsApps\python.exe"
        "$env:LOCALAPPDATA\Microsoft\WindowsApps\python3.exe"
    )

    $removed = 0
    foreach ($stub in $stubs) {
        if (Test-Path $stub) {
            try {
                Remove-Item $stub -Force -ErrorAction Stop
                $removed++
                Write-Host "  Removed Windows Store stub: $stub" -ForegroundColor Green
            } catch {
                Write-Host "  Could not remove $stub - please disable in Settings > Apps > App execution aliases" -ForegroundColor Yellow
            }
        }
    }
    return $removed
}

function Install-PythonViaUv {
    # Use uv to install Python - this puts it in ~/.local/bin which we add to PATH
    $uvPath = "$env:USERPROFILE\.local\bin\uv.exe"
    if (-not (Test-Path $uvPath)) {
        # Try to find uv elsewhere
        $uvPath = (Get-Command uv -ErrorAction SilentlyContinue).Source
    }

    if ($uvPath -and (Test-Path $uvPath)) {
        Write-Host "  Installing Python 3.11 via uv..." -ForegroundColor Yellow
        $result = & $uvPath python install 3.11 2>&1
        if ($LASTEXITCODE -eq 0) {
            Write-Host "  Python 3.11 installed successfully" -ForegroundColor Green
            return $true
        } else {
            Write-Host "  Failed to install Python via uv" -ForegroundColor Red
            return $false
        }
    }
    return $false
}

function Test-Tool {
    param([string]$CheckCommand)
    try {
        $ErrorActionPreference = "SilentlyContinue"
        $output = Invoke-Expression $CheckCommand 2>&1
        $ErrorActionPreference = "Stop"
        if ($LASTEXITCODE -eq 0 -or $output -match "version") {
            return $output | Select-Object -First 1
        }
        return $null
    } catch {
        return $null
    }
}

function Test-IsAdmin {
    $principal = New-Object Security.Principal.WindowsPrincipal([Security.Principal.WindowsIdentity]::GetCurrent())
    return $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
}

function Get-ExecutablePaths {
    param(
        [string]$CommandName,
        [string[]]$FallbackPaths = @()
    )

    $paths = @()
    try {
        $cmd = Get-Command $CommandName -ErrorAction SilentlyContinue | Where-Object { $_.CommandType -eq "Application" } | Select-Object -First 1
        if ($cmd -and $cmd.Source) {
            $paths += $cmd.Source
        }
    } catch {
        # ignore
    }

    foreach ($path in $FallbackPaths) {
        if ($path -and (Test-Path $path)) {
            $paths += $path
        }
    }

    return $paths | Select-Object -Unique
}

function Ensure-WindowsFirewallRules {
    $isWindows = $env:OS -eq "Windows_NT"
    if (-not $isWindows) {
        Write-Host "Windows Firewall rules are only supported on Windows." -ForegroundColor Yellow
        return
    }

    if (-not (Test-IsAdmin)) {
        Write-Host "Windows Firewall updates require an elevated PowerShell." -ForegroundColor Yellow
        Write-Host "Re-run as Administrator: .\\repo.ps1 tools --firewall" -ForegroundColor Cyan
        return
    }

    $homeRoots = @($env:HOME, $env:USERPROFILE) | Where-Object { $_ -and (Test-Path $_) } | Select-Object -Unique
    $minioFallbacks = @()
    foreach ($root in $homeRoots) {
        $minioFallbacks += (Join-Path $root ".sbdev\\bin\\minio.exe")
    }
    $minioFallbacks += (Join-Path $RootDir ".devstack\\relay\\bin\\minio.exe")
    $minioFallbacks += (Join-Path $RootDir "syftbox\\.devstack\\relay\\bin\\minio.exe")

    $nodeFallbacks = @(
        "$env:ProgramFiles\\nodejs\\node.exe",
        "$env:ProgramFiles(x86)\\nodejs\\node.exe"
    )

    $nodePaths = Get-ExecutablePaths -CommandName "node" -FallbackPaths $nodeFallbacks
    $minioPaths = Get-ExecutablePaths -CommandName "minio" -FallbackPaths $minioFallbacks

    if ($nodePaths.Count -eq 0 -and $minioPaths.Count -eq 0) {
        Write-Host "No Node.js or MinIO executables found for firewall rules." -ForegroundColor Yellow
        return
    }

    $ruleCount = 0
    foreach ($path in $nodePaths) {
        $existing = Get-NetFirewallApplicationFilter -Program $path -ErrorAction SilentlyContinue
        if (-not $existing) {
            $name = "BioVault Node.js (" + (Split-Path $path -Leaf) + ")"
            New-NetFirewallRule -DisplayName $name -Direction Inbound -Action Allow -Program $path -Profile Private | Out-Null
            Write-Host "Added firewall rule for Node.js: $path" -ForegroundColor Green
            $ruleCount++
        }
    }

    foreach ($path in $minioPaths) {
        $existing = Get-NetFirewallApplicationFilter -Program $path -ErrorAction SilentlyContinue
        if (-not $existing) {
            $name = "BioVault MinIO (" + (Split-Path $path -Leaf) + ")"
            New-NetFirewallRule -DisplayName $name -Direction Inbound -Action Allow -Program $path -Profile Private | Out-Null
            Write-Host "Added firewall rule for MinIO: $path" -ForegroundColor Green
            $ruleCount++
        }
    }

    if ($ruleCount -eq 0) {
        Write-Host "Firewall rules already present for detected executables." -ForegroundColor Green
    }
}

function Maybe-Configure-WindowsFirewall {
    param(
        [bool]$ConfigureFirewall = $false,
        [bool]$AutoInstall = $false
    )

    if ($ConfigureFirewall) {
        Write-Host ""
        Write-Host "Configuring Windows Firewall rules..." -ForegroundColor Cyan
        Ensure-WindowsFirewallRules
        return
    }

    if ($AutoInstall -and $env:OS -eq "Windows_NT") {
        Write-Host ""
        Write-Host "Add Windows Firewall rules for Node.js and MinIO? [y/N] " -ForegroundColor Cyan -NoNewline
        $response = Read-Host
        if ($response -match "^[Yy]") {
            Ensure-WindowsFirewallRules
        }
    }
}

function Show-Tools {
    param(
        [bool]$AutoInstall = $false,
        [bool]$ConfigureFirewall = $false
    )

    # Check for Windows Store Python stubs that interfere with real Python
    $pythonStubs = @(
        "$env:LOCALAPPDATA\Microsoft\WindowsApps\python.exe"
        "$env:LOCALAPPDATA\Microsoft\WindowsApps\python3.exe"
    )
    $hasStubs = $false
    foreach ($stub in $pythonStubs) {
        if (Test-Path $stub) {
            $hasStubs = $true
            break
        }
    }

    if ($hasStubs) {
        Write-Host "=== Windows Store Python Stubs Detected ===" -ForegroundColor Yellow
        Write-Host "  These stubs intercept Python commands and cause issues." -ForegroundColor Yellow
        Write-Host ""

        if ($AutoInstall) {
            Write-Host "  Removing stubs..." -ForegroundColor Cyan
            Remove-WindowsStorePythonStubs | Out-Null
        } else {
            Write-Host "  Remove them? [Y/n] " -ForegroundColor Cyan -NoNewline
            $response = Read-Host
            if ($response -eq "" -or $response -match "^[Yy]") {
                Remove-WindowsStorePythonStubs | Out-Null
            }
        }
        Write-Host ""
    }

    Write-Host "=== Required Development Tools ===" -ForegroundColor Cyan
    Write-Host ""

    $missingWinget = @()
    $missingCargo = @()

    # Check winget-installed tools
    $versionMismatch = @()
    foreach ($tool in $RequiredTools) {
        $version = Test-Tool $tool.Check
        $reqLabel = if ($tool.Required) { "" } else { " (optional)" }
        if ($version) {
            # Check if version matches required version (if specified)
            $versionOk = $true
            if ($tool.RequiredVersion) {
                if ($version -notmatch $tool.RequiredVersion) {
                    $versionOk = $false
                    $versionMismatch += $tool
                }
            }
            if ($versionOk) {
                Write-Host "  " -NoNewline
                Write-Host "[OK]" -ForegroundColor Green -NoNewline
                Write-Host " $($tool.Name) - $version$reqLabel"
            } else {
                Write-Host "  " -NoNewline
                Write-Host "[!!]" -ForegroundColor Yellow -NoNewline
                Write-Host " $($tool.Name) - $version" -NoNewline
                Write-Host " (need $($tool.RequiredVersion).x)" -ForegroundColor Yellow
            }
        } else {
            # Add all missing tools (both required and optional)
            $missingWinget += $tool
            Write-Host "  " -NoNewline
            Write-Host "[X]" -ForegroundColor Red -NoNewline
            Write-Host " $($tool.Name) - " -NoNewline
            Write-Host "not found" -ForegroundColor Yellow -NoNewline
            Write-Host " ($($tool.Desc))$reqLabel"
        }
    }

    Write-Host ""
    Write-Host "=== Cargo Tools (for pipelines) ===" -ForegroundColor Cyan
    Write-Host ""

    # Check cargo-installed tools
    foreach ($tool in $CargoTools) {
        $version = Test-Tool $tool.Check
        if ($version) {
            Write-Host "  " -NoNewline
            Write-Host "[OK]" -ForegroundColor Green -NoNewline
            Write-Host " $($tool.Name) - $version"
        } else {
            $missingCargo += $tool
            Write-Host "  " -NoNewline
            Write-Host "[X]" -ForegroundColor Red -NoNewline
            Write-Host " $($tool.Name) - " -NoNewline
            Write-Host "not found" -ForegroundColor Yellow -NoNewline
            Write-Host " ($($tool.Desc))"
        }
    }

    Write-Host ""

    $totalMissing = $missingWinget.Count + $missingCargo.Count
    $requiredMissing = ($missingWinget | Where-Object { $_.Required }).Count + $missingCargo.Count

    # Show version mismatch warnings
    if ($versionMismatch.Count -gt 0) {
        Write-Host ""
        Write-Host "=== Version Issues ===" -ForegroundColor Yellow
        foreach ($tool in $versionMismatch) {
            Write-Host "  $($tool.Name): installed version is incompatible" -ForegroundColor Yellow
            if ($tool.ManualInstall) {
                Write-Host "    Download correct version: $($tool.ManualInstall)" -ForegroundColor Cyan
            }
        }
        Write-Host ""
    }

    if ($totalMissing -eq 0 -and $versionMismatch.Count -eq 0) {
        Write-Host "All tools are installed!" -ForegroundColor Green
        Maybe-Configure-WindowsFirewall -ConfigureFirewall $ConfigureFirewall -AutoInstall $AutoInstall
        return
    }

    if ($totalMissing -eq 0 -and $versionMismatch.Count -gt 0) {
        Write-Host "All tools installed but some have version issues (see above)." -ForegroundColor Yellow
        Maybe-Configure-WindowsFirewall -ConfigureFirewall $ConfigureFirewall -AutoInstall $AutoInstall
        return
    }

    if ($requiredMissing -eq 0) {
        Write-Host "All required tools are installed! (some optional tools missing)" -ForegroundColor Green
    }

    # Get unique winget packages (only those with winget defined)
    $wingetPkgs = $missingWinget | Where-Object { $_.Winget -ne $null } | ForEach-Object { $_.Winget } | Sort-Object -Unique

    Write-Host "Missing $totalMissing tool(s)." -ForegroundColor Yellow
    Write-Host ""

    # Show manual install instructions for tools without winget
    $manualTools = $missingWinget | Where-Object { $_.Winget -eq $null -and $_.ManualInstall -ne $null -and $_.UseUv -ne $true }
    if ($manualTools.Count -gt 0) {
        Write-Host "Manual installation required:" -ForegroundColor Yellow
        foreach ($tool in $manualTools) {
            Write-Host "  $($tool.Name): $($tool.ManualInstall)" -ForegroundColor Cyan
        }
        Write-Host ""
    }

    # If --install flag or user confirms, install missing tools
    $doInstall = $AutoInstall
    if (-not $doInstall -and $totalMissing -gt 0) {
        Write-Host "Would you like to install missing tools now? [Y/n] " -ForegroundColor Cyan -NoNewline
        $response = Read-Host
        $doInstall = ($response -eq "" -or $response -match "^[Yy]")
    }

    if ($doInstall) {
        # First, check if any tools need Python via uv
        $needsPython = $missingWinget | Where-Object { $_.UseUv -eq $true }
        if ($needsPython.Count -gt 0) {
            Write-Host ""
            Write-Host "Installing Python via uv..." -ForegroundColor Cyan

            # Make sure uv is available first
            $uvAvailable = Test-Tool "uv --version"
            if (-not $uvAvailable) {
                Write-Host "  -> Installing uv first..." -ForegroundColor Yellow
                $result = cmd /c "winget install --accept-source-agreements --accept-package-agreements -e --id astral-sh.uv 2>&1"
                # Refresh PATH to find uv
                $env:PATH = "$env:USERPROFILE\.local\bin;$env:PATH"
            }

            if (Install-PythonViaUv) {
                # Refresh PATH to find the new Python
                $env:PATH = "$env:USERPROFILE\.local\bin;$env:PATH"
            }
        }

        # Install winget packages
        $wingetTools = $missingWinget | Where-Object { $_.Winget -ne $null }
        if ($wingetTools.Count -gt 0) {
            Write-Host ""
            Write-Host "Installing winget packages..." -ForegroundColor Cyan

            # Track which packages we've already installed (some tools share packages)
            $installedPkgs = @{}

            foreach ($tool in $wingetTools) {
                $pkg = $tool.Winget

                # Skip if we already installed this package
                if ($installedPkgs.ContainsKey($pkg)) {
                    continue
                }
                $installedPkgs[$pkg] = $true

                Write-Host "  -> Installing $pkg..." -ForegroundColor Yellow

                # Check if tool has custom winget args
                if ($tool.WingetArgs) {
                    $installCmd = "winget install --accept-source-agreements --accept-package-agreements -e --id $pkg $($tool.WingetArgs)"
                } else {
                    $installCmd = "winget install --accept-source-agreements --accept-package-agreements -e --id $pkg"
                }

                $result = cmd /c "$installCmd 2>&1"
                if ($LASTEXITCODE -eq 0) {
                    Write-Host "     done" -ForegroundColor Green
                } else {
                    Write-Host "     failed (may need manual install)" -ForegroundColor Red
                    Write-Host "     $result" -ForegroundColor Gray
                }
            }
        }

        # Install cargo packages (only if cargo is available)
        if ($missingCargo.Count -gt 0) {
            $cargoAvailable = Test-Tool "cargo --version"
            if ($cargoAvailable) {
                Write-Host ""
                Write-Host "Installing cargo packages..." -ForegroundColor Cyan
                foreach ($tool in $missingCargo) {
                    Write-Host "  -> Installing $($tool.Crate)..." -ForegroundColor Yellow
                    $result = cmd /c "cargo install $($tool.Crate) --locked 2>&1"
                    if ($LASTEXITCODE -eq 0) {
                        Write-Host "     done" -ForegroundColor Green
                    } else {
                        Write-Host "     failed" -ForegroundColor Red
                    }
                }
            } else {
                Write-Host ""
                Write-Host "Cargo not available - install Rust first, then re-run 'tools --install'" -ForegroundColor Yellow
            }
        }

        Write-Host ""
        Write-Host "Installation complete!" -ForegroundColor Green
        Write-Host ""
        Write-Host "IMPORTANT: Please restart your terminal for PATH changes to take effect." -ForegroundColor Magenta
        Write-Host "Then run '.\repo.ps1 tools' again to verify." -ForegroundColor Magenta
        Maybe-Configure-WindowsFirewall -ConfigureFirewall $ConfigureFirewall -AutoInstall $AutoInstall
    } else {
        # Just show manual install commands
        $uvTools = $missingWinget | Where-Object { $_.UseUv -eq $true }
        $wingetTools = $missingWinget | Where-Object { $_.Winget -ne $null }

        if ($uvTools.Count -gt 0) {
            Write-Host "Python (install via uv):" -ForegroundColor Yellow
            Write-Host "  uv python install 3.11" -ForegroundColor Blue
            Write-Host ""
        }

        if ($wingetTools.Count -gt 0) {
            Write-Host "Winget packages to install:" -ForegroundColor Yellow

            # Track which packages we've shown
            $shownPkgs = @{}

            foreach ($tool in $wingetTools) {
                $pkg = $tool.Winget
                if ($shownPkgs.ContainsKey($pkg)) { continue }
                $shownPkgs[$pkg] = $true

                if ($tool.WingetArgs) {
                    Write-Host "  winget install $pkg $($tool.WingetArgs)" -ForegroundColor Blue
                } else {
                    Write-Host "  winget install $pkg" -ForegroundColor Blue
                }
            }
        }

        if ($missingCargo.Count -gt 0) {
            Write-Host ""
            Write-Host "Cargo packages to install (after Rust is installed):" -ForegroundColor Yellow
            foreach ($tool in $missingCargo) {
                Write-Host "  cargo install $($tool.Crate) --locked" -ForegroundColor Blue
            }
        }

        Write-Host ""
        Write-Host "Or run: .\repo.ps1 tools --install" -ForegroundColor Gray
        if ($ConfigureFirewall) {
            Maybe-Configure-WindowsFirewall -ConfigureFirewall $ConfigureFirewall -AutoInstall $AutoInstall
        }
    }
}

function Get-ManifestProjects {
    if (-not (Test-Path $ManifestPath)) {
        Write-Host "Manifest not found: $ManifestPath" -ForegroundColor Red
        exit 1
    }

    [xml]$manifest = Get-Content $ManifestPath
    $defaultRemote = $manifest.manifest.default.remote
    $defaultRevision = $manifest.manifest.default.revision

    $projects = @()
    foreach ($project in $manifest.manifest.project) {
        $path = if ($project.path) { $project.path } else { $project.name }
        $revision = if ($project.revision) { $project.revision } else { $defaultRevision }
        $remote = if ($project.remote) { $project.remote } else { $defaultRemote }

        # Get remote fetch URL
        $fetchUrl = ""
        foreach ($r in $manifest.manifest.remote) {
            if ($r.name -eq $remote) {
                $fetchUrl = $r.fetch
                break
            }
        }

        $projects += @{
            Name = $project.name
            Path = $path
            Revision = $revision
            Remote = $remote
            FetchUrl = $fetchUrl
        }
    }

    return $projects
}

function Get-HttpsUrl {
    param([string]$SshUrl, [string]$ProjectName)

    # Convert ssh://git@github.com/ to https://github.com/
    if ($SshUrl -match "ssh://git@([^/]+)/") {
        $gitHost = $Matches[1]
        return "https://$gitHost/$ProjectName.git"
    }
    # Convert git@github.com: to https://github.com/
    if ($SshUrl -match "git@([^:]+):") {
        $gitHost = $Matches[1]
        return "https://$gitHost/$ProjectName.git"
    }
    return $SshUrl
}

function Get-SshUrl {
    param([string]$FetchUrl, [string]$ProjectName)

    # Convert to git@host:project.git format
    if ($FetchUrl -match "ssh://git@([^/]+)/") {
        $gitHost = $Matches[1]
        return "git@${gitHost}:$ProjectName.git"
    }
    if ($FetchUrl -match "https?://([^/]+)/") {
        $gitHost = $Matches[1]
        return "git@${gitHost}:$ProjectName.git"
    }
    if ($FetchUrl -match "git@([^:]+):") {
        $gitHost = $Matches[1]
        return "git@${gitHost}:$ProjectName.git"
    }
    return $FetchUrl
}

function Initialize-Workspace {
    param([bool]$UseHttps = $false)

    Write-Host "Initializing workspace..." -ForegroundColor Cyan

    $projects = Get-ManifestProjects

    foreach ($proj in $projects) {
        $targetPath = Join-Path $RootDir $proj.Path

        Write-Host "  -> $($proj.Path): " -NoNewline

        if (Test-Path $targetPath) {
            Write-Host "exists" -ForegroundColor Green
            continue
        }

        # Create parent directory if needed
        $parentDir = Split-Path -Parent $targetPath
        if (-not (Test-Path $parentDir)) {
            New-Item -ItemType Directory -Path $parentDir -Force | Out-Null
        }

        # Determine clone URL
        $cloneUrl = if ($UseHttps) {
            Get-HttpsUrl -SshUrl $proj.FetchUrl -ProjectName $proj.Name
        } else {
            Get-SshUrl -FetchUrl $proj.FetchUrl -ProjectName $proj.Name
        }

        Write-Host "cloning... " -NoNewline -ForegroundColor Yellow

        # Git outputs progress to stderr, so we need to handle this specially
        # Use cmd /c to prevent PowerShell from treating stderr as error
        $null = cmd /c "git clone -b `"$($proj.Revision)`" `"$cloneUrl`" `"$targetPath`" 2>&1"
        if ($LASTEXITCODE -eq 0) {
            Write-Host "done" -ForegroundColor Green
        } else {
            # Try without branch if it fails (might be a SHA)
            $null = cmd /c "git clone `"$cloneUrl`" `"$targetPath`" 2>&1"
            if ($LASTEXITCODE -eq 0) {
                Push-Location $targetPath
                git checkout $proj.Revision 2>&1 | Out-Null
                Pop-Location
                Write-Host "done" -ForegroundColor Green
            } else {
                Write-Host "failed" -ForegroundColor Red
            }
        }
    }

    # Run post-init script if exists
    $postInitPs1 = Join-Path $RootDir "scripts\setup-repo-workspace.ps1"
    $postInitSh = Join-Path $RootDir "scripts\setup-repo-workspace.sh"
    $gitBash = "C:\Program Files\Git\bin\bash.exe"

    if (Test-Path $postInitPs1) {
        Write-Host ""
        Write-Host "Running post-init setup..." -ForegroundColor Cyan
        & $postInitPs1
    } elseif ((Test-Path $postInitSh) -and (Test-Path $gitBash)) {
        Write-Host ""
        Write-Host "Running post-init setup (via bash)..." -ForegroundColor Cyan
        $unixPath = $postInitSh -replace '\\', '/' -replace '^([A-Za-z]):', '/$1'
        $unixPath = $unixPath.ToLower()
        & $gitBash $unixPath
    }

    Write-Host ""
    Write-Host "Workspace initialized!" -ForegroundColor Green
}

function Sync-Workspace {
    Write-Host "Syncing workspace..." -ForegroundColor Cyan

    $projects = Get-ManifestProjects

    foreach ($proj in $projects) {
        $targetPath = Join-Path $RootDir $proj.Path

        Write-Host "  -> $($proj.Path): " -NoNewline

        if (-not (Test-Path $targetPath)) {
            Write-Host "missing (run --init first)" -ForegroundColor Red
            continue
        }

        Push-Location $targetPath
        try {
            git fetch --all 2>&1 | Out-Null
            $ErrorActionPreference = "SilentlyContinue"
            $output = git checkout $proj.Revision 2>&1
            $ErrorActionPreference = "Stop"
            if ($LASTEXITCODE -eq 0) {
                Write-Host "synced to $($proj.Revision)" -ForegroundColor Green
            } else {
                Write-Host "failed" -ForegroundColor Red
            }
        } finally {
            Pop-Location
        }
    }
}

function Get-RepoBranch {
    param([string]$RepoPath)
    Push-Location $RepoPath
    try {
        $ErrorActionPreference = "SilentlyContinue"
        $branch = git rev-parse --abbrev-ref HEAD 2>$null
        if ($LASTEXITCODE -eq 0 -and $branch) {
            return $branch.Trim()
        }
        return "HEAD"
    } finally {
        $ErrorActionPreference = "Stop"
        Pop-Location
    }
}

function Test-RepoDirty {
    param([string]$RepoPath)
    Push-Location $RepoPath
    try {
        $ErrorActionPreference = "SilentlyContinue"
        $status = git status --porcelain -uno 2>$null
        return ($status -and $status.Length -gt 0)
    } finally {
        $ErrorActionPreference = "Stop"
        Pop-Location
    }
}

function Test-RepoDirtyAll {
    param([string]$RepoPath)
    Push-Location $RepoPath
    try {
        $ErrorActionPreference = "SilentlyContinue"
        $status = git status --porcelain 2>$null
        return ($status -and $status.Length -gt 0)
    } finally {
        $ErrorActionPreference = "Stop"
        Pop-Location
    }
}

function Get-RepoShortSha {
    param([string]$RepoPath)
    Push-Location $RepoPath
    try {
        $ErrorActionPreference = "SilentlyContinue"
        $sha = git rev-parse --short HEAD 2>$null
        if ($sha) {
            return $sha.Trim()
        }
        return "unknown"
    } finally {
        $ErrorActionPreference = "Stop"
        Pop-Location
    }
}

function Show-Tree {
    # Show root repo status
    $rootBranch = Get-RepoBranch $RootDir
    $rootDirty = if (Test-RepoDirty $RootDir) { " [dirty]" } else { "" }
    $rootName = Split-Path -Leaf $RootDir

    Write-Host "$rootName/ " -ForegroundColor Cyan -NoNewline
    Write-Host "[$rootBranch]" -ForegroundColor Blue -NoNewline
    if ($rootDirty) {
        Write-Host $rootDirty -ForegroundColor Red
    } else {
        Write-Host ""
    }

    $projects = Get-ManifestProjects

    foreach ($proj in $projects) {
        $targetPath = Join-Path $RootDir $proj.Path
        $name = Split-Path -Leaf $proj.Path

        Write-Host "  +-- " -NoNewline
        Write-Host "$name/ " -ForegroundColor Cyan -NoNewline

        if (-not (Test-Path $targetPath)) {
            Write-Host "[missing]" -ForegroundColor Red
            continue
        }

        $branch = Get-RepoBranch $targetPath
        $dirty = Test-RepoDirty $targetPath
        $manifestRev = $proj.Revision

        # Show branch or detached state
        if ($branch -eq "HEAD") {
            $sha = Get-RepoShortSha $targetPath
            Write-Host "(detached: $sha)" -ForegroundColor Yellow -NoNewline
        } else {
            Write-Host "[$branch]" -ForegroundColor Blue -NoNewline
        }

        # Show dirty status
        if ($dirty) {
            Write-Host " [dirty]" -ForegroundColor Red -NoNewline
        }

        # Show if different from manifest
        if ($branch -ne "HEAD" -and $branch -ne $manifestRev) {
            Write-Host " [manifest:$manifestRev]" -ForegroundColor Blue -NoNewline
        }

        Write-Host ""
    }

    Write-Host ""
    Write-Host "Legend:" -ForegroundColor Green
    Write-Host "  [branch]     - on branch" -ForegroundColor Blue
    Write-Host "  (detached)   - detached HEAD" -ForegroundColor Yellow
    Write-Host "  [dirty]      - uncommitted changes" -ForegroundColor Red
    Write-Host "  [manifest:x] - differs from manifest revision" -ForegroundColor Blue
    Write-Host "  [missing]    - repo not checked out" -ForegroundColor Red
}

function Invoke-FetchAll {
    Write-Host "Fetching all repos..." -ForegroundColor Cyan

    $projects = Get-ManifestProjects

    foreach ($proj in $projects) {
        $targetPath = Join-Path $RootDir $proj.Path

        Write-Host "  -> $($proj.Path): " -NoNewline

        if (-not (Test-Path $targetPath)) {
            Write-Host "missing" -ForegroundColor Red
            continue
        }

        Push-Location $targetPath
        try {
            $ErrorActionPreference = "SilentlyContinue"
            $output = git fetch --all --prune 2>&1
            $ErrorActionPreference = "Stop"
            if ($LASTEXITCODE -eq 0) {
                $branch = Get-RepoBranch $targetPath
                if ($branch -eq "HEAD") {
                    Write-Host "fetched (detached)" -ForegroundColor Yellow
                } else {
                    # Check ahead/behind
                    $ErrorActionPreference = "SilentlyContinue"
                    $upstream = git rev-parse --abbrev-ref "@{u}" 2>$null
                    $ErrorActionPreference = "Stop"
                    if ($LASTEXITCODE -eq 0 -and $upstream) {
                        $ErrorActionPreference = "SilentlyContinue"
                        $counts = git rev-list --left-right --count "HEAD...$upstream" 2>$null
                        $ErrorActionPreference = "Stop"
                        if ($counts -match "(\d+)\s+(\d+)") {
                            $ahead = [int]$Matches[1]
                            $behind = [int]$Matches[2]
                            if ($ahead -eq 0 -and $behind -eq 0) {
                                Write-Host "up-to-date" -ForegroundColor Green
                            } else {
                                $status = @()
                                if ($ahead -gt 0) { $status += "ahead:$ahead" }
                                if ($behind -gt 0) { $status += "behind:$behind" }
                                Write-Host ($status -join " ") -ForegroundColor Yellow
                            }
                        } else {
                            Write-Host "fetched" -ForegroundColor Green
                        }
                    } else {
                        Write-Host "fetched (no upstream)" -ForegroundColor Yellow
                    }
                }
            } else {
                Write-Host "failed" -ForegroundColor Red
            }
        } finally {
            Pop-Location
        }
    }
}

function Invoke-PullAll {
    param([bool]$UseRebase = $false)

    Write-Host "Pulling all repos..." -ForegroundColor Cyan

    $projects = Get-ManifestProjects

    foreach ($proj in $projects) {
        $targetPath = Join-Path $RootDir $proj.Path

        Write-Host "  -> $($proj.Path): " -NoNewline

        if (-not (Test-Path $targetPath)) {
            Write-Host "missing" -ForegroundColor Red
            continue
        }

        Push-Location $targetPath
        try {
            $branch = Get-RepoBranch $targetPath

            if ($branch -eq "HEAD") {
                Write-Host "detached (skipped)" -ForegroundColor Yellow
                continue
            }

            if (Test-RepoDirty $targetPath) {
                Write-Host "dirty (skipped)" -ForegroundColor Yellow
                continue
            }

            $ErrorActionPreference = "SilentlyContinue"
            if ($UseRebase) {
                $output = git pull --rebase 2>&1
            } else {
                $output = git pull 2>&1
            }
            $ErrorActionPreference = "Stop"

            if ($LASTEXITCODE -eq 0) {
                Write-Host "updated" -ForegroundColor Green
            } else {
                Write-Host "failed" -ForegroundColor Red
            }
        } finally {
            Pop-Location
        }
    }
}

function Invoke-Lint {
    param([string[]]$Args)

    $force = $false
    $verbose = $false

    foreach ($arg in $Args) {
        switch ($arg) {
            "--force" { $force = $true }
            "--verbose" { $verbose = $true }
            "--help" { Write-Host "Usage: .\\repo.ps1 lint [--force] [--verbose]" -ForegroundColor Yellow; return }
            "-h" { Write-Host "Usage: .\\repo.ps1 lint [--force] [--verbose]" -ForegroundColor Yellow; return }
            default {
                Write-Host "Unknown lint option: $arg" -ForegroundColor Red
                return
            }
        }
    }

    $startTime = Get-Date
    $dirtyRepos = @()
    $projects = Get-ManifestProjects

    if ($force) {
        if (Test-Path (Join-Path $RootDir "lint.sh")) {
            $dirtyRepos += "."
        }
        foreach ($proj in $projects) {
            $repoPath = Join-Path $RootDir $proj.Path
            if (Test-Path (Join-Path $repoPath ".git")) {
                $dirtyRepos += $proj.Path
            }
        }
    } else {
        if (Test-RepoDirtyAll $RootDir) {
            $dirtyRepos += "."
        }
        foreach ($proj in $projects) {
            $repoPath = Join-Path $RootDir $proj.Path
            if (Test-Path $repoPath) {
                if (Test-RepoDirtyAll $repoPath) {
                    $dirtyRepos += $proj.Path
                }
            }
        }

        if ($dirtyRepos.Count -eq 0) {
            Write-Host "No dirty repos" -ForegroundColor Green
            return
        }
    }

    $tmpDir = Join-Path ([System.IO.Path]::GetTempPath()) ("repo-lint-" + [System.Guid]::NewGuid().ToString("N"))
    New-Item -ItemType Directory -Path $tmpDir | Out-Null

    $winPs1Path = Join-Path $RootDir "win.ps1"
    if (-not (Test-Path $winPs1Path)) {
        Write-Host "win.ps1 not found at $winPs1Path" -ForegroundColor Red
        return
    }

    $jobs = @()
    $reposWithLint = @()

    try {
        foreach ($repo in $dirtyRepos) {
            $repoPath = if ($repo -eq ".") { $RootDir } else { Join-Path $RootDir $repo }
            $label = if ($repo -eq ".") { "biovault-desktop" } else { $repo }
            $repoKey = if ($repo -eq ".") { "root" } else { $repo }
            $scriptPath = Join-Path $repoPath "lint.sh"

            if (-not (Test-Path $scriptPath)) {
                continue
            }

            $reposWithLint += $repo
            Write-Host "-> $label" -ForegroundColor Cyan

            $outfile = Join-Path $tmpDir (($repoKey -replace '[\\/]', '_') + ".out")
            $errfile = Join-Path $tmpDir (($repoKey -replace '[\\/]', '_') + ".err")
            $start = Get-Date

            $job = Start-Job -ScriptBlock {
                param($jobRepoPath, $jobWinPs1Path, $jobOutFile, $jobErrFile)
                Push-Location $jobRepoPath
                try {
                    & $jobWinPs1Path ./lint.sh 1> $jobOutFile 2> $jobErrFile
                    return $LASTEXITCODE
                } finally {
                    Pop-Location
                }
            } -ArgumentList $repoPath, $winPs1Path, $outfile, $errfile

            $jobs += [pscustomobject]@{
                Repo = $repo
                Label = $label
                Job = $job
                OutFile = $outfile
                ErrFile = $errfile
                Start = $start
            }
        }

        if ($reposWithLint.Count -eq 0) {
            Write-Host "No repos with lint.sh" -ForegroundColor Green
            return
        }

        Write-Host ""
        $failures = 0
        $passed = 0

        foreach ($item in $jobs) {
            Wait-Job -Job $item.Job | Out-Null

            $elapsed = [int]((Get-Date) - $item.Start).TotalSeconds
            $elapsedText = " (${elapsed}s)"

            $jobResult = Receive-Job -Job $item.Job -ErrorAction SilentlyContinue
            $exitCode = if ($jobResult -is [System.Array]) { $jobResult[-1] } else { $jobResult }
            if ($null -eq $exitCode) { $exitCode = 1 }

            if ($exitCode -ne 0) {
                Write-Host "FAIL $($item.Label)$elapsedText" -ForegroundColor Red
                if (Test-Path $item.OutFile) { Get-Content $item.OutFile }
                if (Test-Path $item.ErrFile) { Get-Content $item.ErrFile }
                Write-Host ""
                $failures++
            } elseif ($verbose) {
                Write-Host "OK $($item.Label)$elapsedText" -ForegroundColor Green
                if (Test-Path $item.OutFile) { Get-Content $item.OutFile }
                if (Test-Path $item.ErrFile) { Get-Content $item.ErrFile }
                Write-Host ""
                $passed++
            } else {
                Write-Host "OK $($item.Label)$elapsedText" -ForegroundColor Green
                $passed++
            }

            Remove-Job -Job $item.Job -Force | Out-Null
        }

        $totalTime = [int]((Get-Date) - $startTime).TotalSeconds
        Write-Host ""
        if ($failures -eq 0) {
            Write-Host "$passed repos passed in ${totalTime}s" -ForegroundColor Green
        } else {
            Write-Host "$failures failed, $passed passed in ${totalTime}s" -ForegroundColor Red
            exit 1
        }
    } finally {
        if (Test-Path $tmpDir) {
            Remove-Item -Recurse -Force $tmpDir
        }
    }
}

function Invoke-Test {
    param([string[]]$Args)

    $winPs1Path = Join-Path $RootDir "win.ps1"
    if (-not (Test-Path $winPs1Path)) {
        Write-Host "win.ps1 not found at $winPs1Path" -ForegroundColor Red
        return
    }

    & $winPs1Path ./repo test @Args
    if ($LASTEXITCODE -ne 0) {
        exit $LASTEXITCODE
    }
}

function Invoke-SshRewrite {
    Write-Host "Rewriting remotes to SSH..." -ForegroundColor Cyan

    $projects = Get-ManifestProjects

    foreach ($proj in $projects) {
        $targetPath = Join-Path $RootDir $proj.Path

        Write-Host "  -> $($proj.Path): " -NoNewline

        if (-not (Test-Path $targetPath)) {
            Write-Host "missing" -ForegroundColor Red
            continue
        }

        Push-Location $targetPath
        try {
            $ErrorActionPreference = "SilentlyContinue"
            $currentUrl = git remote get-url origin 2>$null
            $ErrorActionPreference = "Stop"

            if ($currentUrl -match "^git@" -or $currentUrl -match "^ssh://") {
                Write-Host "ok (already SSH)" -ForegroundColor Green
                continue
            }

            $sshUrl = Get-SshUrl -FetchUrl $proj.FetchUrl -ProjectName $proj.Name
            $ErrorActionPreference = "SilentlyContinue"
            git remote set-url origin $sshUrl 2>$null
            $ErrorActionPreference = "Stop"

            if ($LASTEXITCODE -eq 0) {
                Write-Host "updated" -ForegroundColor Green
            } else {
                Write-Host "failed" -ForegroundColor Red
            }
        } finally {
            Pop-Location
        }
    }
}

function Invoke-CheckoutMain {
    Write-Host "Checking out main in all repos..." -ForegroundColor Cyan
    Invoke-Checkout -Rev "main" -Target "all" -DoReset $false
}

function Resolve-Targets {
    param([string[]]$Targets)

    $projects = Get-ManifestProjects
    $resolved = @()

    foreach ($target in $Targets) {
        if ($target -eq "all") {
            $resolved += $projects | ForEach-Object { $_.Path }
        } elseif ($target -eq "self") {
            $resolved += "."
        } else {
            # Find by path or name
            $found = $projects | Where-Object {
                $_.Path -eq $target -or
                (Split-Path -Leaf $_.Path) -eq $target -or
                $_.Name -eq $target
            }
            if ($found) {
                $resolved += $found.Path
            } else {
                Write-Host "Unknown target: $target" -ForegroundColor Red
            }
        }
    }

    return $resolved | Sort-Object -Unique
}

function Invoke-Checkout {
    param(
        [string]$Rev,
        [string]$Target,
        [bool]$DoReset = $false
    )

    if (-not $Rev -or -not $Target) {
        Write-Host "Usage: .\repo.ps1 checkout [--reset] <rev> <target>" -ForegroundColor Red
        return
    }

    $targets = Resolve-Targets @($Target)

    foreach ($path in $targets) {
        $targetPath = if ($path -eq ".") { $RootDir } else { Join-Path $RootDir $path }
        $displayPath = if ($path -eq ".") { "self" } else { $path }

        Write-Host "  -> ${displayPath}: " -NoNewline

        if (-not (Test-Path $targetPath)) {
            Write-Host "missing" -ForegroundColor Red
            continue
        }

        Push-Location $targetPath
        try {
            if ($DoReset) {
                git reset --hard 2>&1 | Out-Null
                git clean -fd 2>&1 | Out-Null
            }

            git fetch --all 2>&1 | Out-Null
            $ErrorActionPreference = "SilentlyContinue"
            $output = git checkout $Rev 2>&1
            $ErrorActionPreference = "Stop"

            if ($LASTEXITCODE -eq 0) {
                Write-Host "checked out" -ForegroundColor Green
            } else {
                Write-Host "failed" -ForegroundColor Red
            }
        } finally {
            Pop-Location
        }
    }
}

function Invoke-Switch {
    param(
        [string]$Branch,
        [string[]]$Targets,
        [bool]$Create = $false
    )

    if (-not $Branch) {
        Write-Host "Usage: .\repo.ps1 switch [-b] <branch> <targets...>" -ForegroundColor Red
        return
    }

    if ($Targets.Count -eq 0) {
        Write-Host "Missing targets. Use 'all' or list repo names/paths." -ForegroundColor Red
        return
    }

    $resolvedTargets = Resolve-Targets $Targets

    foreach ($path in $resolvedTargets) {
        $targetPath = if ($path -eq ".") { $RootDir } else { Join-Path $RootDir $path }
        $displayPath = if ($path -eq ".") { "self" } else { $path }

        Write-Host "  -> ${displayPath}: " -NoNewline

        if (-not (Test-Path $targetPath)) {
            Write-Host "missing" -ForegroundColor Red
            continue
        }

        Push-Location $targetPath
        try {
            if (Test-RepoDirty $targetPath) {
                Write-Host "dirty (skipped)" -ForegroundColor Yellow
                continue
            }

            # Check if branch exists locally
            $ErrorActionPreference = "SilentlyContinue"
            git show-ref --verify --quiet "refs/heads/$Branch" 2>$null
            $localExists = ($LASTEXITCODE -eq 0)
            $ErrorActionPreference = "Stop"

            if ($localExists) {
                $ErrorActionPreference = "SilentlyContinue"
                $output = git checkout $Branch 2>&1
                $ErrorActionPreference = "Stop"
                if ($LASTEXITCODE -eq 0) {
                    Write-Host "switched" -ForegroundColor Green
                } else {
                    Write-Host "failed" -ForegroundColor Red
                }
                continue
            }

            # Check if branch exists on remote
            $ErrorActionPreference = "SilentlyContinue"
            $remoteBranch = git ls-remote --heads origin $Branch 2>$null
            $ErrorActionPreference = "Stop"
            if ($remoteBranch) {
                $ErrorActionPreference = "SilentlyContinue"
                $output = git checkout -b $Branch "origin/$Branch" 2>&1
                $ErrorActionPreference = "Stop"
                if ($LASTEXITCODE -eq 0) {
                    Write-Host "switched (from remote)" -ForegroundColor Green
                } else {
                    Write-Host "failed" -ForegroundColor Red
                }
                continue
            }

            # Create new branch if -b flag
            if ($Create) {
                $ErrorActionPreference = "SilentlyContinue"
                $output = git checkout -b $Branch 2>&1
                $ErrorActionPreference = "Stop"
                if ($LASTEXITCODE -eq 0) {
                    Write-Host "created" -ForegroundColor Green
                } else {
                    Write-Host "failed" -ForegroundColor Red
                }
                continue
            }

            Write-Host "branch not found" -ForegroundColor Red
        } finally {
            Pop-Location
        }
    }
}

function Update-ManifestRevision {
    param(
        [string]$Target,
        [string]$Revision
    )

    if (-not $Target -or -not $Revision) {
        Write-Host "Usage: .\repo.ps1 track <target> <rev>" -ForegroundColor Red
        return
    }

    $projects = Get-ManifestProjects
    $matches = @($projects | Where-Object {
        $_.Path -eq $Target -or
        (Split-Path -Leaf $_.Path) -eq $Target
    })

    if ($matches.Count -eq 0) {
        Write-Host "Project not found: $Target" -ForegroundColor Red
        return
    }

    if ($matches.Count -gt 1) {
        Write-Host "Ambiguous target. Use full path." -ForegroundColor Red      
        return
    }

    $targetPath = $matches[0].Path
    $targetName = $matches[0].Name

    $manifestText = Get-Content $ManifestPath -Raw
    $updated = $false

    $escapedPath = [regex]::Escape($targetPath)
    $patternPath = '(<project\\b[^>]*\\bpath="' + $escapedPath + '"[^>]*\\brevision=")([^"]*)(")'
    if ([regex]::IsMatch($manifestText, $patternPath)) {
        $manifestText = [regex]::Replace(
            $manifestText,
            $patternPath,
            { param($m) "$($m.Groups[1].Value)$Revision$($m.Groups[3].Value)" },
            1
        )
        $updated = $true
    }

    if (-not $updated) {
        $escapedName = [regex]::Escape($targetName)
        $patternName = '(<project\\b[^>]*\\bname="' + $escapedName + '"[^>]*\\brevision=")([^"]*)(")'
        if ([regex]::IsMatch($manifestText, $patternName)) {
            $manifestText = [regex]::Replace(
                $manifestText,
                $patternName,
                { param($m) "$($m.Groups[1].Value)$Revision$($m.Groups[3].Value)" },
                1
            )
            $updated = $true
        }
    }

    if (-not $updated) {
        Write-Host "Project not found in manifest: $targetPath" -ForegroundColor Red
        return
    }

    $utf8NoBom = New-Object System.Text.UTF8Encoding($false)
    [System.IO.File]::WriteAllText($ManifestPath, $manifestText, $utf8NoBom)
    Write-Host "Updated $ManifestFile - $targetPath -> $Revision" -ForegroundColor Green
}

# Parse command line arguments manually (to match bash --option style)
$Command = ""
$UseHttps = $false
$UseRebase = $false
$DoReset = $false
$CreateBranch = $false
$AutoInstall = $false
$ConfigureFirewall = $false
$PositionalArgs = @()

for ($i = 0; $i -lt $args.Count; $i++) {
    $arg = $args[$i]
    switch -Regex ($arg) {
        "^(--help|-h)$" {
            Show-Usage
            exit 0
        }
        "^--init$" {
            $Command = "init"
            break
        }
        "^--https$" {
            $UseHttps = $true
            break
        }
        "^--rebase$" {
            $UseRebase = $true
            break
        }
        "^--reset$" {
            $DoReset = $true
            break
        }
        "^--install$" {
            $AutoInstall = $true
            break
        }
        "^--firewall$" {
            $ConfigureFirewall = $true
            break
        }
        "^--force$" {
            if ($Command -eq "lint" -or $Command -eq "test") {
                $PositionalArgs += $arg
                break
            }
            Write-Host "Unknown option: $arg" -ForegroundColor Red
            Show-Usage
            exit 1
        }
        "^--verbose$" {
            if ($Command -eq "lint") {
                $PositionalArgs += $arg
                break
            }
            Write-Host "Unknown option: $arg" -ForegroundColor Red
            Show-Usage
            exit 1
        }
        "^--help$" {
            if ($Command -eq "lint" -or $Command -eq "test") {
                $PositionalArgs += $arg
                break
            }
            Write-Host "Unknown option: $arg" -ForegroundColor Red
            Show-Usage
            exit 1
        }
        "^-b$" {
            $CreateBranch = $true
            break
        }
        "^--.+" {
            # Unknown --option
            Write-Host "Unknown option: $arg" -ForegroundColor Red
            Show-Usage
            exit 1
        }
        default {
            if (-not $Command) {
                $Command = $arg
            } else {
                $PositionalArgs += $arg
            }
        }
    }
}

# Route to appropriate command
switch ($Command) {
    "" {
        Show-Tree
    }
    "init" {
        Initialize-Workspace -UseHttps $UseHttps
    }
    "sync" {
        Sync-Workspace
    }
    "tools" {
        Show-Tools -AutoInstall $AutoInstall -ConfigureFirewall $ConfigureFirewall
    }
    "fetch" {
        Invoke-FetchAll
    }
    "pull" {
        Invoke-PullAll -UseRebase $UseRebase
    }
    "lint" {
        Invoke-Lint -Args $PositionalArgs
    }
    "test" {
        Invoke-Test -Args $PositionalArgs
    }
    "ssh" {
        Invoke-SshRewrite
    }
    "main" {
        Invoke-CheckoutMain
    }
    "checkout" {
        if ($PositionalArgs.Count -ge 2) {
            Invoke-Checkout -Rev $PositionalArgs[0] -Target $PositionalArgs[1] -DoReset $DoReset
        } else {
            Write-Host "Usage: .\repo.ps1 checkout [--reset] <rev> <target>" -ForegroundColor Red
        }
    }
    "switch" {
        if ($PositionalArgs.Count -ge 2) {
            $branch = $PositionalArgs[0]
            $targets = $PositionalArgs[1..($PositionalArgs.Count - 1)]
            Invoke-Switch -Branch $branch -Targets $targets -Create $CreateBranch
        } else {
            Write-Host "Usage: .\repo.ps1 switch [-b] <branch> <targets...>" -ForegroundColor Red
        }
    }
    "track" {
        if ($PositionalArgs.Count -ge 2) {
            Update-ManifestRevision -Target $PositionalArgs[0] -Revision $PositionalArgs[1]
        } else {
            Write-Host "Usage: .\repo.ps1 track <target> <rev>" -ForegroundColor Red
        }
    }
    default {
        # Check if it's a two-arg shortcut for track (e.g., "repo.ps1 biovault main")
        if ($PositionalArgs.Count -ge 1) {
            Update-ManifestRevision -Target $Command -Revision $PositionalArgs[0]
        } else {
            Write-Host "Unknown command: $Command" -ForegroundColor Red
            Show-Usage
            exit 1
        }
    }
}
