# setup-podman-windows.ps1 - Setup Podman (Hyper-V) on Windows for BioVault
# Usage: .\scripts\setup-podman-windows.ps1

$ErrorActionPreference = "Stop"

Write-Host "=== System Info ===" -ForegroundColor Cyan
Write-Host "OS: $([System.Environment]::OSVersion.VersionString)"
whoami

Write-Host ""
Write-Host "=== Hyper-V Status ===" -ForegroundColor Cyan
try {
    Get-WindowsOptionalFeature -Online -FeatureName Microsoft-Hyper-V* | Format-Table FeatureName, State
} catch {
    Write-Host "Could not query Hyper-V features. If this fails, run PowerShell as Administrator." -ForegroundColor Yellow
}

Write-Host ""
Write-Host "=== Installing Chocolatey (if needed) ===" -ForegroundColor Cyan
if (-not (Get-Command choco -ErrorAction SilentlyContinue)) {
    Set-ExecutionPolicy Bypass -Scope Process -Force
    [System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072
    Invoke-Expression ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1'))
    $env:PATH = "$env:PATH;C:\ProgramData\chocolatey\bin"
}
choco --version

Write-Host ""
Write-Host "=== Installing Podman CLI (if needed) ===" -ForegroundColor Cyan
if (-not (Get-Command podman -ErrorAction SilentlyContinue)) {
    choco install podman-cli -y
}
$env:PATH = "C:\ProgramData\chocolatey\bin;$env:PATH"

Write-Host ""
Write-Host "=== Podman Version ===" -ForegroundColor Cyan
podman --version

Write-Host ""
Write-Host "=== Initializing Podman Machine ===" -ForegroundColor Cyan

# Compute resource targets
$totalMemMb = [int](([float](Get-ComputerInfo).OsTotalVisibleMemorySize) / 1024)
$desiredMemMb = 16384
if ($totalMemMb - 512 -lt $desiredMemMb) {
    $desiredMemMb = [Math]::Max(2048, $totalMemMb - 512)
    Write-Host "Requested 16GB exceeds available memory. Using ${desiredMemMb}MB instead." -ForegroundColor Yellow
}
$logicalCpus = (Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfLogicalProcessors -Sum).Sum
$desiredCpus = 8
if ($logicalCpus -lt $desiredCpus) {
    $desiredCpus = [Math]::Max(1, [int]$logicalCpus)
    Write-Host "Requested 8 CPUs exceeds available logical CPUs. Using ${desiredCpus} CPU(s) instead." -ForegroundColor Yellow
}

function Ensure-PodmanMachine {
    param(
        [string]$Provider,
        [string]$Name
    )

    $env:CONTAINERS_MACHINE_PROVIDER = $Provider
    [Environment]::SetEnvironmentVariable("CONTAINERS_MACHINE_PROVIDER", $Provider, "Process")

    # Keep existing machine if already created
    $existingMachines = @()
    try { $existingMachines = podman machine list --format json | ConvertFrom-Json } catch { $existingMachines = @() }
    $hasMachine = $false
    foreach ($m in $existingMachines) {
        if ($m.Name -eq $Name) { $hasMachine = $true }
    }
    if (-not $hasMachine) {
        podman machine init $Name --cpus $desiredCpus --memory $desiredMemMb
        if ($LASTEXITCODE -ne 0) { throw "podman machine init failed for $Name ($Provider)" }
    }

    podman machine start $Name
    if ($LASTEXITCODE -ne 0) {
        $state = ""
        try { $state = podman machine inspect $Name --format '{{.State}}' } catch { $state = "" }
        if ($state -ne "running") { throw "podman machine start failed for $Name ($Provider)" }
    }
    podman system connection default $Name
    if ($LASTEXITCODE -ne 0) { throw "podman system connection default failed for $Name ($Provider)" }
}

try {
    Write-Host "Attempting Hyper-V Podman machine..." -ForegroundColor Cyan
    Ensure-PodmanMachine -Provider "hyperv" -Name "podman-hyperv"
} catch {
    Write-Host "Hyper-V setup failed; falling back to WSL provider." -ForegroundColor Yellow
    Ensure-PodmanMachine -Provider "wsl" -Name "podman-wsl"

    $wslConfigPath = "$env:USERPROFILE\.wslconfig"
    if (-not (Test-Path $wslConfigPath)) {
        $memSetting = "${desiredMemMb}MB"
        $wslConfig = @"
[wsl2]
memory=$memSetting
processors=$desiredCpus
"@
        Set-Content -Path $wslConfigPath -Value $wslConfig -Encoding UTF8
        Write-Host "Created $wslConfigPath to set WSL resources (memory=$memSetting, processors=$desiredCpus)." -ForegroundColor Yellow
        Write-Host "Run 'wsl --shutdown' to apply WSL resource changes." -ForegroundColor Yellow
    } else {
        Write-Host "WSL config already exists at $wslConfigPath; please ensure it sets memory/processors as desired." -ForegroundColor Yellow
    }
}

Write-Host ""
Write-Host "=== Testing Podman ===" -ForegroundColor Cyan
podman info
podman run --rm hello-world

Write-Host ""
Write-Host "Podman setup complete. You can now run: .\\win.ps1 ./test-scenario.sh --docker tests/scenarios/syqure-distributed.yaml" -ForegroundColor Green
