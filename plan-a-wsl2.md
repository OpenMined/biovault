## Plan A: WSL2 on Windows CI (Namespace)

Goal: prove WSL2 works on the runner quickly, then enable the HERC2 nested
pipeline using Podman (WSL backend) with minimal CI iteration.

### Fast Local Baseline (already working here)
1) Confirm WSL2 + Podman local path
   - `wsl --status`
   - `wsl -l -v` (verify version 2)
   - `podman machine init` (WSL provider auto)
   - `podman machine start`
2) Run quick local tests
   - `./win.ps1 tests/scripts/quick-nextflow-test.sh`
   - `./win.ps1 tests/scripts/quick-herc2-test.sh`

### Minimal CI Probe Workflow (independent GH action)
Purpose: keep the probe under 1-2 minutes to iterate fast.

1) Print version/status
   - `wsl --status`
   - `wsl --list --verbose`
2) Try WSL2 install / upgrade
   - `wsl --set-default-version 2`
   - `wsl --install -d Ubuntu-24.04` (or Debian)
3) Verify distro runs
   - `wsl -d Ubuntu-24.04 --exec uname -a`
   - `wsl -d Ubuntu-24.04 --exec cat /etc/os-release`

### Fallback: WSL2 without Store
If `wsl --install` fails on CI, use `wsl --import` with a rootfs tar:
1) Download rootfs (Ubuntu or Alpine)
2) `wsl --import <DistroName> C:\wsl\<DistroName> <rootfs.tar>`
3) `wsl --set-version <DistroName> 2`
4) `wsl -d <DistroName> --exec uname -a`

### Optional Debugging on CI
- Use `nsc rdp <instance-id>` to inspect runner state interactively.
- Use a reverse SSH tunnel (tunshell or equivalent) only if needed for quick
  iterative debugging. Keep it short-lived and minimal.

### If WSL2 Probe Succeeds
1) Add a second CI job that runs:
   - `./win.ps1 tests/scripts/quick-nextflow-test.sh`
2) If that passes, run the HERC2 quick test:
   - `./win.ps1 tests/scripts/quick-herc2-test.sh`
3) Only then expand to the full scenario tests.
