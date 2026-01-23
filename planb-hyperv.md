## Plan B: Hyper-V Podman VM Copy Mode (No Host Mounts)

Goal: run HERC2 nested containers on Windows CI using Hyper-V even if host
mounts fail due to 9P/junction issues.

### Phase 0: Confirm Hyper-V Provider Locally
1) Force Hyper-V
   - `set CONTAINERS_MACHINE_PROVIDER=hyperv` (Git Bash: `export ...`)
2) Verify provider
   - `podman machine inspect --format "{{.VMType}}"`
3) Start the VM
   - `podman machine start`

### Phase 1: Test "Clean Path" Mounts (No Junctions)
Hypothesis: mounts work if the host path has no junctions.

1) Move test inputs to a junction-free path, e.g. `C:\bv-data`
2) Check the VM can see it
   - `podman machine ssh -- ls /mnt/c/bv-data`
3) Run quick tests
   - `./win.ps1 tests/scripts/quick-nextflow-test.sh`
   - `./win.ps1 tests/scripts/test-hyperv-herc2.sh`

If this works, we can avoid VM copy mode and just standardize on clean paths
for CI inputs.

### Phase 2: VM Copy Mode (No Host Mounts)
If mounts still fail, copy the project + data into the VM and run everything
from VM-local paths.

#### Copy In (Git Bash example)
1) Package inputs + project
   - `tar -C <project_root> -cf - . | podman machine ssh -- mkdir -p /tmp/bv-project && tar -xf - -C /tmp/bv-project`
   - `tar -C <data_root> -cf - . | podman machine ssh -- mkdir -p /tmp/bv-data && tar -xf - -C /tmp/bv-data`

#### Run Inside the VM
2) Execute Nextflow in the VM (no host mounts)
   - `podman machine ssh -- bash -lc "cd /tmp/bv-project && nextflow run workflow.nf -with-podman -work-dir /tmp/bv-work"`
   - If using `bv`, call the CLI inside the VM only if it is installed there.
     Otherwise, run Nextflow directly.

#### Copy Out Results
3) Pull results back to Windows
   - `podman machine ssh -- tar -C /tmp/bv-results -cf - . | tar -C <host_results_dir> -xf -`

Notes:
- Use `/tmp` in the VM for Nextflow temp/logs to avoid 9P/SMB issues.
- Keep the test dataset small for fast iteration.

### Phase 3: Automate for CI
1) Wrap the copy-run-copy sequence in a script:
   - `tests/scripts/hyperv-vm-copy-herc2.sh`
2) Use the 1-row samplesheet for quick feedback in CI.
3) Only scale up once the small run is stable.
