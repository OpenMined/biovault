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
Try Hyper-V approach:
```powershell
# Enable Hyper-V (may need runner support)
Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Hyper-V -All -NoRestart

# Install Podman
choco install podman-cli -y

# Set provider to Hyper-V
$env:CONTAINERS_MACHINE_PROVIDER = "hyperv"

# Initialize and start
podman machine init
podman machine start

# Create docker alias
Set-Alias -Name docker -Value podman
```

### Fallback
If Hyper-V doesn't work on CI:
- Run Docker-dependent tests (`import-and-run.yaml`) on Linux runners only
- Run non-Docker tests (`messaging-core.yaml`, `key-management.yaml`) on Windows

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

## Next Steps

1. ~~**Test Hyper-V on CI**: Enable Hyper-V feature and try Podman with hyperv provider~~ ✅ Done!
2. ~~**Nested container support**: Enable Nextflow to spawn task containers with Podman~~ ✅ Done!
3. **Verify full pipeline execution**: Test `import-and-run.yaml` scenario with Podman
4. **Document for users**: Update user docs for Podman as Docker alternative

## References

- [Podman for Windows](https://github.com/containers/podman/blob/main/docs/tutorials/podman-for-windows.md)
- [Vampire/setup-wsl Action](https://github.com/Vampire/setup-wsl)
- [WSL Installation](https://docs.microsoft.com/en-us/windows/wsl/install)
