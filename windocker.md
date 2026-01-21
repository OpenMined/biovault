# Running Linux Containers on Windows CI

## Problem Statement
The `import-and-run.yaml` scenario requires Nextflow, which pulls Linux container images. Windows Server Docker Engine only supports Windows containers natively - it cannot run Linux containers without additional setup.

## Approaches Tried

### 1. Docker Desktop with DockerCli.exe Switch
**Status**: ❌ Not available on CI runners

```powershell
$dockerCli = "C:\Program Files\Docker\Docker\DockerCli.exe"
& $dockerCli -SwitchLinuxEngine
```

- Docker Desktop isn't installed on namespace runners
- Only Windows Docker Engine is available
- `docker info` shows `OSType: windows`

### 2. WSL2 with Docker Inside
**Status**: ❌ WSL not enabled on namespace runners

```powershell
wsl --status  # Returns "Access is denied"
wsl --install -d Debian  # Error: 0x8000ffff
```

- Namespace runners run as `nt authority\system` with `BUILTIN\Administrators`
- WSL feature is not enabled at the runner level
- `wsl --set-default-version 2` works (just sets registry)
- Actually using WSL distributions fails with access denied

### 3. Vampire/setup-wsl GitHub Action
**Status**: ❌ Permission denied

```yaml
- uses: Vampire/setup-wsl@v6
  with:
    distribution: Ubuntu-24.04
    wsl-version: 2
```

Error:
```
WslRegisterDistribution failed with error: 0x80070005
Error: 0x80070005 Access is denied.
```

### 4. Podman with WSL Provider
**Status**: ❌ WSL not available

```powershell
choco install podman-cli -y
podman machine init  # Auto-detects WSL
```

Error:
```
Error: command C:\Windows\system32\wsl.exe [-l --quiet] failed: exit status 0xffffffff ()
```

- `podman-cli` from Chocolatey works for installation
- But `podman machine init` fails because WSL isn't available
- The `--provider` flag isn't available in `podman-cli` (remote client only)

### 5. Podman with Hyper-V Provider
**Status**: ✅ Working on CI!

```powershell
$env:CONTAINERS_MACHINE_PROVIDER = "hyperv"
podman machine init
podman machine start
```

- Hyper-V is already enabled on namespace runners
- No reboot required
- Podman 5.7.1 works with Hyper-V backend
- Successfully runs Linux containers (Fedora CoreOS 43 VM)
- Docker socket compatibility via `/var/run/docker.sock` symlink

### 6. Rancher Desktop
**Status**: ❌ Didn't work on CI

- Tried installing via Chocolatey
- GUI application, hard to run headless
- Requires WSL or Hyper-V backend anyway

## Technical Details

### Namespace Runner Environment
```
OS: Microsoft Windows NT 10.0.20348.0
User: nt authority\system
Groups: BUILTIN\Administrators (enabled)
Runner: namespace-profile-windows-medium
```

### Windows Docker Engine Info
```
OSType: windows
Storage Driver: windowsfilter
Default Isolation: process
```

### What Works
- Windows Docker Engine for Windows containers
- `docker pull` with `--platform=linux/amd64` flag (but can't run)
- Chocolatey package installation
- PowerShell with admin privileges

### What Doesn't Work
- WSL commands (access denied at feature level)
- Docker Linux containers without WSL or Hyper-V
- Registering WSL distributions
- Direct Hyper-V mounts of junctioned paths (use host-mount mode)

## Podman Configuration

### Provider Selection
Podman provider can be set via:
1. **Environment variable**: `CONTAINERS_MACHINE_PROVIDER=hyperv` or `wsl`
2. **Config file**: `%APPDATA%\containers\containers.conf`
   ```toml
   [machine]
   provider = "hyperv"
   ```
3. **MSI installer**: Select during installation

### Podman-cli vs Full Podman
- `podman-cli` (Chocolatey): Remote client only, no `--provider` flag
- Full Podman MSI: Has machine management with provider selection
- Both can use environment variable for provider

## Recommended Solution

### For CI (namespace runners without WSL)
Use Podman Hyper-V with host-mount mode:

```powershell
# Hyper-V is already enabled on namespace runners
choco install podman-cli -y

$env:CONTAINERS_MACHINE_PROVIDER = "hyperv"
$env:BIOVAULT_CONTAINER_RUNTIME = "podman"
$env:BIOVAULT_HYPERV_MOUNT = "1"
$env:BIOVAULT_HYPERV_HOST_DIR = "C:vtemp"  # optional

podman machine init
podman machine start

Set-Alias -Name docker -Value podman
```

## Local Development (Windows with WSL2)

On a Windows machine with WSL2 available:
```powershell
# Install via Chocolatey
choco install podman-cli -y

# Initialize (auto-detects WSL)
podman machine init
podman machine start

# Test
podman run --rm hello-world
```

This works because WSL2 is properly enabled on desktop Windows.

## Files Modified

- `biovault/.github/workflows/biovault-scenario-tests.yml` - CI workflow
- `biovault/.github/workflows/ci.yml` - Main CI entry point
- `biovault/cli/src/cli/commands/run_dynamic.rs` - Container runtime auto-detection (Docker/Podman)
- `biovault/cli/src/config.rs` - Added `container_runtime` config option

## Container Runtime Configuration

The CLI now supports both Docker and Podman with automatic detection:

### Auto-detection (default)
1. Tries configured Docker path (from `binary_paths.docker`)
2. Tries system `docker` command
3. Falls back to `podman` if Docker isn't available or isn't running Linux containers

### Environment Variable
```bash
BIOVAULT_CONTAINER_RUNTIME=podman  # Force Podman
BIOVAULT_CONTAINER_RUNTIME=docker  # Force Docker
BIOVAULT_CONTAINER_RUNTIME=auto    # Auto-detect (default)
```

### Config File
In `config.yaml`:
```yaml
binary_paths:
  container_runtime: podman  # or "docker" or "auto"
```

## Nested Container Support (Nextflow Tasks)

When running Nextflow pipelines, the outer Nextflow container needs to spawn inner task containers. This requires:

### Docker Desktop Mode
- Mounts `/var/run/docker.sock` into the Nextflow container
- Uses `docker.enabled = true` in Nextflow config
- Task containers use Docker CLI to talk to Docker daemon

### Podman Mode
- Mounts Podman socket (e.g., `/run/user/1000/podman/podman.sock`) into container
- Sets `CONTAINER_HOST=unix:///run/podman/podman.sock` environment variable
- Uses `podman.enabled = true` in Nextflow config
- Uses `/bin/sh` shell for alpine-based containers (no bash)
- The `nextflow-runner` image includes both Docker CLI and Podman CLI

### Nextflow Runner Image

The `ghcr.io/openmined/nextflow-runner:25.10.2` image contains:
- Nextflow 25.10.2
- Docker CLI 28.0.1 (for Docker Desktop mode)
- Podman CLI 5.3.1 (for Podman mode)

Build from `docker/windows/Dockerfile.nextflow-runner`:
```dockerfile
FROM nextflow/nextflow:25.10.2
# ... installs Docker CLI and Podman CLI
```

### Runtime Config Generation

The CLI automatically generates a `.biovault-runtime.config` file in the project directory:

**For Docker:**
```groovy
process.executor = 'local'
docker.enabled = true
docker.runOptions = '-u $(id -u):$(id -g)'
```

**For Podman:**
```groovy
process.executor = 'local'
podman.enabled = true
process.shell = ['/bin/sh', '-ue']
```

This config is passed to Nextflow with the `-c` flag.

## Shell Compatibility Issues

### The `printf '%q'` Problem

The bioscript workflow templates used `printf '%q'` for shell-safe quoting:
```groovy
script:
def genoFileName = genotype_file.getName()
"""
GENO_FILE=\$(printf '%q' "${genoFileName}")
bioscript classify "${assets_dir}/classify_herc2.py" --file \$GENO_FILE --participant_id "${participant_id}"
"""
```

This fails with `invalid directive` error because:
- Nextflow containers use `/bin/sh` (configured via `process.shell = ['/bin/sh', '-ue']`)
- `printf '%q'` is a **bash-specific** extension, not available in POSIX `/bin/sh`
- Alpine-based containers (like bioscript) only have `/bin/sh`

### The Fix

Remove the `printf '%q'` wrapper and quote variables directly:
```groovy
script:
def genoFileName = genotype_file.getName()
"""
bioscript classify "${assets_dir}/classify_herc2.py" --file "${genoFileName}" --participant_id "${participant_id}"
"""
```

This is safe because:
- Nextflow's `def genoFileName = genotype_file.getName()` already gives a clean filename
- The filename comes from pipeline config, not untrusted user input
- Double-quoting in shell handles spaces and special characters

### Affected Files

The fix was applied to all workflow files in the bioscript repo:
- `examples/herc2/herc2-classifier/workflow.nf` (already on main)
- `examples/apol1/apol1-classifier/workflow.nf`
- `examples/brca/brca-classifier/workflow.nf`
- `examples/thalassemia/thalassemia-classifier/workflow.nf`
- `python/src/bioscript/biovault.py` (the workflow generator)

PR: `fix/shell-compatibility` branch in bioscript repo

## Path Format Handling

### CSV Path Requirements

BioVault CLI extracts file paths from CSV/TSV inputs to set up container mounts. For proper mount extraction:

| Path Format | Works With | Example |
|-------------|------------|---------|
| Windows (`C:/...`) | BioVault CLI mount detection | `C:/Users/admin/data/file.vcf` |
| WSL (`/mnt/c/...`) | Podman internal paths | `/mnt/c/Users/admin/data/file.vcf` |
| Docker (`/c/...`) | Docker Desktop internal paths | `/c/Users/admin/data/file.vcf` |

**Important**: Use Windows-style paths (`C:/Users/...`) in CSV files so the CLI can properly detect and mount the directories.

### Path Conversion in Test Scripts

The quick test script (`tests/scripts/quick-herc2-test.sh`) includes a helper:
```bash
# Convert Git Bash path (/c/Users/...) to Windows path (C:/Users/...)
windows_path() {
    local p="$1"
    echo "$p" | sed 's|^/\([a-zA-Z]\)/|\1:/|'
}
```

## Hyper-V vs WSL2 Filesystem Issues

### The Problem

Podman can use two backends on Windows:
- **WSL2**: Uses 9P/Plan9 protocol for Windows mounts, better POSIX compatibility
- **Hyper-V**: Uses different mount mechanism (SMB/CIFS-like), stricter limitations

On CI runners (namespace-profile-windows), WSL2 is not available, so we use Hyper-V. This causes I/O errors when Nextflow tries to perform POSIX operations on Windows-mounted paths:

```
mv: preserving times for '/mnt/c/.../nxf-tmp.XXX': Invalid argument
mv: listing attributes of '/mnt/c/.../nxf-tmp.XXX': Invalid argument
ERROR ~ Input/output error
```

### The Fix

Set environment variables to redirect Nextflow temp files to native Linux paths:

```rust
// In run_dynamic.rs - Podman-specific settings
.arg("-e").arg("NXF_HOME=/tmp/.nextflow")  // Nextflow home directory
.arg("-e").arg("NXF_TEMP=/tmp")             // Temp files (nxf-tmp.*)
```

Also use `/tmp` for the Nextflow log file:
```rust
let docker_log_path = if using_podman {
    "/tmp/.nextflow.log".to_string()
} else {
    windows_path_to_container(&nextflow_log_path, using_podman)
};
```

### Why It Works Locally But Fails on CI

| Environment | Backend | Windows Mount Protocol | POSIX Compatibility |
|-------------|---------|------------------------|---------------------|
| Local dev machine | WSL2 | 9P/Plan9 | Good - most operations work |
| CI runners | Hyper-V | 9P (Fedora CoreOS VM) | **Broken** - permission denied |

The fix ensures all Nextflow internal files go to `/tmp` (native Linux filesystem inside the Podman VM), while actual workflow data still uses mounted paths.

### Hyper-V 9P Mount Issue (Workaround)

The Hyper-V backend uses a 9P file sharing stack that struggles with Windows
directory junctions (common under user profiles). Directly mounting repo paths
can produce "Permission denied" or empty directories inside the VM, especially
for nested containers.

Current workaround: use Hyper-V host-mount mode to stage inputs into a
junction-free directory (e.g., `C:vtemp`) and rewrite CSV paths to that flat
root. This avoids junctions and keeps nested containers working.

### Hyper-V Host Mount Mode (Junction-Free)

To work around Hyper-V 9P and junction issues, the CLI supports a host-mount
mode that stages inputs into a junction-free directory and rewrites CSV paths.
This keeps all mounted paths under a clean root like `%SystemDrive%\\bvtemp`.

Enable via env vars:
```
BIOVAULT_HYPERV_MOUNT=1
BIOVAULT_HYPERV_HOST_DIR=C:\bvtemp   # optional; default is %SystemDrive%\bvtemp
BIOVAULT_KEEP_HYPERV_HOST_DIR=1      # optional; keep staging dir for debugging
```

In this mode:
- Inputs are copied into a flat `data/` directory under the host mount root.
- CSV paths are rewritten to point at that flat directory.
- Nextflow uses `process.stageInMode = 'copy'` to avoid symlink issues.
- Results are written to a mounted `results/` dir and copied back to the requested output path.

## Local Testing

### Quick Test Script

For rapid local testing without full CI:
```bash
./win.ps1 tests/scripts/quick-herc2-test.sh
```

This script:
1. Sets up required environment variables (`SYC_VAULT`)
2. Creates a temporary samplesheet with Windows paths
3. Runs `biovault run` with the herc2 pipeline
4. Cleans up temporary files

### Successful Test Results

Local testing with Podman nested containers:
```
[c0/3b0bc9] Submitted process > USER:herc2_classifier (p001)
[23/7c95e6] Submitted process > USER:herc2_classifier (p002)
...
[c0/3b0bc9] Completed process > USER:herc2_classifier (p001)
[23/7c95e6] Completed process > USER:herc2_classifier (p002)
```

All 13 herc2 classifier tasks complete successfully with the shell compatibility fix.

## Next Steps

1. Merge bioscript PR: Merge `fix/shell-compatibility` branch to update all workflows
2. Document for users: Update user docs for Podman as Docker alternative
3. Re-enable nightly schedule once CI is stable

## Summary

**Local Development (Windows with WSL2)**:
- OK: Install Podman CLI via Chocolatey
- OK: Use WSL2 backend (default): `podman machine init && podman machine start`
- OK: Nested containers work

**CI (namespace-profile-windows)**:
- WSL not available (permission denied)
- OK: Podman Hyper-V with host-mount mode runs nested containers
- Use `BIOVAULT_HYPERV_MOUNT=1` and a flat host dir (default `%SystemDrive%vtemp`)

## References

- [Podman for Windows](https://github.com/containers/podman/blob/main/docs/tutorials/podman-for-windows.md)
- [Vampire/setup-wsl Action](https://github.com/Vampire/setup-wsl)
- [WSL Installation](https://docs.microsoft.com/en-us/windows/wsl/install)

