# SyftBox Rust Library Migration Plan (BioVault)

## Where BioVault calls the SyftBox *binary* today

- `syftbox-sdk/src/syftbox/control.rs`: spawns `syftbox` via `Command::new(...)`, waits for start/stop, `ps`/`kill` PID management; passes Go-oriented flags like `--control-plane` plus `--client-url/--client-token`.
- `cli/src/cli/commands/daemon.rs`: BioVault daemon calls `start_syftbox()`/`is_syftbox_running()` from `syftbox-sdk` to ensure the external client daemon is running.
- `cli/src/cli/commands/check.rs`: probes `syftbox -v` (and has special sbenv lookup) to report dependency status/version.
- `sbenv/cli/src/main.rs`: manages external SyftBox daemons by spawning binaries (`nohup`, `pgrep`, `kill`, `curl` health checks).
- Scripts/infra that assume a `syftbox` executable exists: `install.sh`, `cli/src/deps.yaml`, `tests/os/test_*_install.*`, `justfile`, `test-scenario.sh` (devstack chooses go|rust client binaries).

## What BioVault actually *needs* from SyftBox at runtime

- BioVault’s “app logic” is filesystem-based (`datasites/<email>/app_data/...` + ACL files) via `syftbox-sdk` helpers; it does not fundamentally depend on the SyftBox daemon being a separate process.
- The only hard coupling to the external binary is: “keep the datasite tree synced” + “optionally provide a local control-plane HTTP port”.

## Target architecture

- Keep `syftbox-sdk` as BioVault’s “datasite layout + RPC file conventions + crypto helpers” crate.
- Replace the *process-spawn* backend in `syftbox-sdk::syftbox::control` with an *embedded* backend that runs the SyftBox Rust daemon (from `syftbox/rust`) as in-process tasks.
- Preserve a fallback “external process” backend for a while (and for `sbenv`/dev tooling), controlled by a flag/env/config.

## Phased port plan

### 1) Turn `syftbox/rust` into a usable library (upstream within the submodule)

- Add a library entrypoint in the package (e.g. `syftbox/rust/src/lib.rs`) and make the daemon pieces public:
  - Export `config::Config`, `control::ControlPlane`, `client::Client`, `workspace`, etc.
  - Move the daemon wiring currently in `syftbox/rust/src/main.rs` into a library function/struct (e.g. `daemon::run(cfg, opts)` or `DaemonHandle::start(...)`).
- Refactor for embedding:
  - Remove/avoid installing a `ctrl_c` handler inside library code (`client::Client::start` currently spawns one); accept an explicit shutdown signal/token so BioVault can own lifecycle.
  - Provide a graceful shutdown path for the control-plane server (axum `with_graceful_shutdown`) and the sync loops.
  - Allow specifying log path (BioVault currently “sandboxes” SyftBox by setting `HOME`/`SYFTBOX_DATA_DIR`; embedded mode should write logs under `${data_dir}/.syftbox/logs/` rather than the user home default).
- Keep `syftbox/rust/src/main.rs` as a thin CLI wrapper that calls the library API (so the binary still works for devstack).

### 2) Add an embedded backend to `syftbox-sdk`

- In `syftbox-sdk`, introduce a backend abstraction in `syftbox-sdk/src/syftbox/control.rs` (or split into `control_process.rs` + `control_embedded.rs`):
  - `ProcessBackend`: current behavior (spawn/ps/kill).
  - `EmbeddedBackend`: starts/stops the in-process daemon via the new `syftbox/rust` library API.
- Add a runtime selector:
  - Env: `BV_SYFTBOX_BACKEND=process|embedded` (or config-driven via BioVault config.yaml later).
  - Default initially remains `process` to reduce risk; flip later.
- Update “is running” detection to be compatible with embedded:
  - Prefer probing the configured `client_url` (TCP connect or `/v1/status`) rather than `ps` matching `syftbox`, because embedded mode won’t show up as a `syftbox` process. This also makes process mode more robust.

### 3) Wire BioVault daemon to use embedded SyftBox

- In `cli/src/cli/commands/daemon.rs`:
  - Start embedded SyftBox once during daemon startup, store a handle inside the BioVault daemon struct, and shut it down on daemon exit.
  - Keep the existing “ensure running / restart attempts” logic, but make it operate on the embedded handle (restart if the task ends) rather than trying to re-spawn an external process.
- Ensure config compatibility:
  - Continue using BioVault’s `SyftboxRuntimeConfig` (`cli/src/config.rs:to_syftbox_runtime_config`) to set `config_path` + `data_dir` + `email`.
  - Map that into the syftbox-rust `ConfigOverrides` so embedded uses BioVault’s data-dir layout (BioVault default is `${BIOVAULT_HOME}` not `~/SyftBox`).

### 4) Deprecate the external `syftbox` prerequisite in BioVault UX (after parity confidence)

- `cli/src/cli/commands/check.rs`: when backend is `embedded`, report SyftBox as “provided by BioVault (embedded)” and show the embedded library version; keep the old `syftbox -v` probe only for `process`.
- `cli/src/cli/commands/setup.rs` / `cli/src/deps.yaml`: stop installing SyftBox when embedded is selected; keep optional install instructions for users who still want the standalone client.
- Update dev/test tooling:
  - Extend `test-scenario.sh` client modes to include `embedded` (or map `rust` → embedded when running `bv daemon`), and update devstack scripts to not assume a `syftbox` executable for that mode.
  - Keep `sbenv` as an external-process manager (it can remain “process-only” by design).

## Key behavior decision to make up front

- Today, BioVault starts SyftBox but does not reliably stop it (because it’s a separate daemon). Embedded mode will tie SyftBox lifetime to the BioVault daemon process; confirm that’s desirable for Desktop (or we add a small dedicated `bv syftboxd` background service that hosts the embedded library in its own process).

