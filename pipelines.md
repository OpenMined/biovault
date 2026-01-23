# BioVault pipelines and projects: deep dive

This document describes how pipelines and projects are represented in BioVault, how they are executed, how Nextflow templates are wired, and how multi-party sharing works. It also covers how the app imports/runs projects and closes with two rewrite recommendations focused on flexible multiparty workflows and reusable sub-units.

## Mental model

BioVault treats a **project** as a self-contained unit of execution with a declared input/output contract, and a **pipeline** as a higher-level DAG that wires multiple projects together. Multiparty execution is handled at the pipeline runner level: a given datasite executes only the steps that target it, while shared outputs are published into SyftBox shared folders with `syft.pub.yaml` permissions.

Core entities:
- **Datasite**: a SyftBox identity + local data root (e.g., `SYFTBOX_EMAIL` and `SYFTBOX_DATA_DIR`).
- **Project**: a versioned workflow with inputs/outputs, typically `dynamic-nextflow` or `shell`.
- **Pipeline**: a list of project steps with explicit bindings between step outputs and step inputs.
- **Submission**: a shared copy of a project (and its assets) in `shared/biovault/submissions/...` plus `syft.pub.yaml` and a project message.
- **Message**: the inbox payload that tells a recipient where a submission lives and how to run it.

## Project specification and execution

### Project files
A project is a folder containing:
- `project.yaml` (spec)
- `workflow.nf` (Nextflow workflow, if `dynamic-nextflow`) # an entry point, which could also be a shell script if the template is shell or something else, so this gives us flexibiltiy to upgrade in the future
- `assets/` (scripts/data shared with the project)

Something missing here is s a version string for these files (and as schema name) like in kubernetes so we can easily version and add features

### Project spec (`project.yaml`)
Key fields:
- `name`, `author`, `version`
- `template`: `dynamic-nextflow` or `shell`
- `workflow`: `workflow.nf` or `workflow.sh`
- `assets`: list of asset file paths to include in submissions
- `inputs` and `outputs`
- `parameters` (user-configurable runtime settings)
- `datasites` (optional recipient list for submissions)

### Input types
The type system is defined in `biovault/cli/src/project_spec.rs` and is shared between project inputs and pipeline inputs. Supported primitives and composites:
- `String`, `Bool`, `File`, `Directory`
- `ParticipantSheet`, `GenotypeRecord`, `BiovaultContext`
- `List[T]`, `Map[String, T]`, `Record{field: Type}`
- `?` suffix for optional (e.g., `File?`)

Type compatibility is validated in the pipeline runner before execution. The pipeline runner checks that step bindings match the expected input types and can resolve literals like `File(path)`.

### Dynamic Nextflow template (`template: dynamic-nextflow`)

Runtime path (as implemented in `biovault/cli/src/cli/commands/run_dynamic.rs`):
1. `project.yaml` is loaded and validated.
2. A `template.nf` is loaded from `~/.biovault/env/dynamic-nextflow/template.nf`.
3. Inputs and parameters are converted into JSON payloads (`inputs.json` and `params.json`).
4. BioVault injects `assets_dir` and `results_dir` into params.
5. Nextflow is launched with the template + workflow.

Importantly here we are using command line arguments and env variables to be able to call any external tool at a step in the pipeline, in this case its nextflow, in other cases it could be a syqure binary or a bash or docker command etc.

User workflow contracts:
- The workflow is named `USER`.
- The first argument to `take:` is always `context` (the `BiovaultContext`).
- The remaining `take:` parameters correspond to `inputs` from `project.yaml`.
- Outputs are emitted by name, matching `project.yaml` outputs.

Type mapping (from `docs/yaml-types-reference.md`):
- `File` -> Nextflow `path`
- `Directory`, `String`, `Bool` -> Nextflow `val`
- Optional `File?` can be empty

The generated context provides:
- `context.params` (parameters + injected `assets_dir` and `results_dir`)
- `context.inputs` (input metadata)

### Shell template (`template: shell`)
Shell projects run a `workflow.sh` script. The runner injects environment variables, including:
- `BV_PROJECT_DIR`, `BV_RESULTS_DIR`, `BV_ASSETS_DIR`
- `BV_INPUT_<INPUT_NAME>` and `BV_OUTPUT_<OUTPUT_NAME>`
- `BV_DATASITES`, `BV_CURRENT_DATASITE`, `BV_DATASITE_INDEX`
- `BV_SYFTBOX_DATA_DIR`, `BV_DATASITES_ROOT`
- `BV_BIN` (path to `bv` if set by the caller)

Shell support only got added today so dont take any of this stuff as expected i think its just highlighting the kinds of metadata that processes potentially need downstream to do their work flexibilty

The shell runner resolves templates inside input/output paths using:
- `{current_datasite}`, `{datasites.index}`, `{datasite.index}`, `{datasites}`

## Pipeline specification and execution

### Pipeline spec (`pipeline.yaml`)
Core fields (see `biovault/cli/src/pipeline_spec.rs` and `docs/pipeline-system-guide.md`):
- `name`
- `inputs` (optional pipeline-level inputs, with optional defaults)
- `steps` (ordered list of project runs)

Step fields:
- `id`: step identifier
- `uses`: project path or registered project name
- `with`: input bindings
- `publish`: output aliases
- `store`: database storage (SQL)
- `runs_on`: datasites to run on (new)
- `foreach`: legacy datasite list (still supported)
- `order`: optional (currently `parallel` just warns and runs sequentially)
- `share`: publish outputs into shared datasite storage (new)

### Step bindings
Bindings support:
- `inputs.<name>` -> pipeline input
- `step.<id>.outputs.<name>` -> upstream step output
- `step.<id>.outputs.<name>.manifest` -> manifest file containing per-datasite outputs
- Literals: `File(path)`, `Directory(path)`, `String(value)`, etc.

When a step output is referenced across multiple datasites, the pipeline runner writes a manifest file under:
```
<results>/manifests/<step>/<output>_paths.txt
```
Each line is `datasite<TAB>path`.

### Multi-datasite execution model
The pipeline runner (in `biovault/cli/src/cli/commands/pipeline.rs`) executes one datasite at a time:
- It resolves the **current datasite** from `BIOVAULT_DATASITE_OVERRIDE`, config email, `SYFTBOX_EMAIL`, or `BIOVAULT_DATASITE`.
- If a step has `runs_on`/`foreach`, only the matching datasite executes it.
- To force a single process to run all targets, set `BIOVAULT_PIPELINE_RUN_ALL=1`.

- i wouldnt use env variables for all this stuff up front i think we can just define this in the yaml spec first and leave the implementation details to the supported runtime / template (the purpose of the template if you check the nextflow ones is to provide execution context and guard rails to running other peoples code by controlling what goes in and where it can run)

For each step run, the runner sets:
- `BIOVAULT_DATASITE_OVERRIDE` to the target datasite.
- `BIOVAULT_DATASITES_OVERRIDE` to the full list of step targets.

`BIOVAULT_DATASITES_OVERRIDE` is used by the shell runner to render `{datasites.index}` and set `BV_DATASITES` consistently even when the project’s own `datasites` list is empty.

### Sharing outputs
The pipeline `share` block defines file sharing at the pipeline level so projects don’t need to manually write `syft.pub.yaml`:

```yaml
share:
  allele_freq_shared:
    source: allele_freq
    path: shared/biovault/shares/{run_id}/{current_datasite}/allele_freq.tsv
    read: [client1@sandbox.local, client2@sandbox.local]
    write: [client1@sandbox.local, client2@sandbox.local]
    admin: [aggregator@sandbox.local]
```

This is an important part of the system which is, we need the built in concept of sharing files between different users and then waiting on them, the code that will copy the file and handle the syft.pub.yaml permission files as well as the downstream code to wait on it should be baked in so these things can be expressed really simply.

Behavior:
- `path` can be a `syft://` URL or a path under the current datasite root.
- The runner writes `syft.pub.yaml` in the parent directory of `path`.
- The shared output is recorded as a `syft://...` URL in step outputs.
- When downstream steps bind to that output, the pipeline runner resolves the `syft://` to a local filesystem path based on `SYFTBOX_DATA_DIR`.

Available template variables for `share.path`:
- `{current_datasite}`
- `{datasites.index}` / `{datasite.index}`
- `{datasites}`
- `{run_id}` (pipeline run ID)

## Submission and inbox flow (app import behavior)

### Submission (`bv submit`)
Submitting a project does all of the following (see `biovault/cli/src/cli/commands/submit.rs`):
1. Copies `project.yaml`, `workflow` and `assets` into a shared folder under the sender datasite.
2. Encrypts assets using SyftBox storage and the recipient list (`project.datasites` or the target recipient).
3. Writes `syft.pub.yaml` with read/write permissions for recipients.
4. Sends a project message containing:
   - `project_location` (a `syft://...` URL)
   - metadata (assets, participants, human text, sender/receiver paths)

syft://{datasite}/path urls are important because they circumvent the need to have DNS/IP of services that are not online while retaining iddntity, and they prevent path traversal hacks by normalizing everything to a datasites/ root and make it easy tot hink about where you are reading and writing from in a networked context rather than a local file system implementation detail.

### Inbox processing (`bv message process`)
Processing a project message (see `biovault/cli/src/cli/commands/messages.rs`) does:
1. Resolves the `syft://` URL to a local path under `SYFTBOX_DATA_DIR`.
2. Copies the submitted project into a local run directory.
3. Executes the project (dynamic-nextflow or shell) in a `results-test` or `results-real` dir.
4. Copies results back into the submission folder and optionally approves (shares results).

This is the primary mechanism by which the app imports and runs projects. Pipelines can be invoked indirectly via a shell project that runs `bv run <pipeline.yaml>`.

## Variables and runtime context summary

### Datasite resolution
Current datasite is resolved using:
1. the yaml has a template to self reference who you are in the computation
2. `Config::load().email`
3. `SYFTBOX_EMAIL`
4. `BIOVAULT_DATASITE`


### Shell project environment
- `BV_INPUT_<NAME>`, `BV_OUTPUT_<NAME>` # inputs
- `BV_PROJECT_DIR`, `BV_RESULTS_DIR`, `BV_ASSETS_DIR` # places where things ar eexeucintg
- `BV_DATASITES`, `BV_CURRENT_DATASITE`, `BV_DATASITE_INDEX` # whos executing them now
- `BV_SYFTBOX_DATA_DIR`, `BV_DATASITES_ROOT`, `BV_BIN` # where things are on ont he system

## Example: share-kitchen-sink (multi-party pipeline)

This scenario exercises:
- A step that runs on all datasites and shares a file.
- A step that runs on client datasites only and shares a derived file.
- An aggregator step that consumes a manifest of shared files.
- A rebroadcast step that shares a combined output back to clients.

See:
- `tests/scenarios/share-kitchen-sink/assets/pipeline.yaml`
- `tests/scenarios/share-kitchen-sink/assets/share-*/project.yaml`

## Notes on reusability

Current reuse model:
- **Local**: `bv project import` registers a project by name; pipelines can `uses: <name>`.
- **Path-based**: `uses: ./relative/path`.
- **Submissions**: `bv submit` shares a project to another datasite (with assets and permissions).

Missing pieces that affect reuse today:
- No built-in version pinning in `uses:`.
- No standard lockfile to record exact commit or hash of the project spec.
- Sharing is local to a datasite (via SyftBox), not a global registry.

## Recommendations if we rewrote this from scratch

1) Introduce a versioned registry + lockfile for reusable sub-units.
   - Allow `uses: ://project@1.2.3` and `uses: git+https://...@sha`. # need to think about how to reference these
   - Store a resolved lockfile in pipelines to freeze input/output schemas, runner type, and asset digests.
   - Make the pipeline runner resolve and cache projects before execution, so users can mix local and remote steps without modifying upstream specs.

2) Make data sharing and multiparty routing first-class in the execution DAG.
   - Model `share`, `collect`, and `broadcast` as built-in step types with typed inputs/outputs and implicit SyftBox permissions.
   - Push these into the pipeline engine rather than shell scripts, so projects can stay single-party and reusable.
   - Expose a minimal runner interface (`run(inputs) -> outputs`) and let runtimes be pluggable (Nextflow, shell, Python, container) with uniform metadata and data movement semantics.

These two changes keep multiparty logic in the pipeline layer, while allowing projects to stay small, reusable, and ETL-focused without requiring downstream edits.


okay so i want to create some good yaml specifications including schema names and versions that cover all this stuff so we can do an upgrade first to make sure things work and a re-usable.

i think we can introduce hashes for reusable compeonents so that they can be resolved locally (as well as making them less strict if you just want to use the dev versio lcaolly)

i want native support for referecnign datasites, doing parallell work in different lcoations or just doing it on one location or some subset, and the ability to do round robin by say referencing next_datasite or something

Heres a reduce operation prototypei made a long time ago that shows how you could model that
author: "madhava@openmined.org"
project: "add"
language: "python"
description: "Add two numbers"
code:
  - functions.py

# Define shared resources using anchors
shared_inputs:
  data: &data FilePipe("{datasite}/data/data.txt")
  output: &output FilePipe("{datasite}/fedreduce/{project}/data/{step}/result.txt")

shared_outputs:
  result: &result FilePipe("{author}/fedreduce/{project}/data/result/result.txt")

# Define the main workflow parameters
workflow:
  datasites: &datasites []

# Define what each step does
steps:
  - first:
      inputs:
        - a: StaticPipe(0) # Override input for the first step
  - last:
      output:
        path: *result
        permissions:
          read:
            - *datasites
  - foreach: *datasites
    run: "{datasite}"
    function: "add"
    inputs:
      - a: FilePipe("{prev_datasite}/fedreduce/{project}/data/{prev_step}/result.txt")
      - b: *data
    output:
      path: *output
      permissions:
        read:
          - "{next_datasite}"

complete:
  exists: *result


we need start and stop conditions, waiting all that stuff.
we need somethings that can be static and others that are depenednent ont he local runtime environent.

I want your best thinking and expertise and lets start with a good initial spec that doesnt back us into a corner.

The idea that you can build a small ETL between steps for your own uses means things could be really re-usable.

First question:
- what should we call the outer unit and what should the composable units be called

Proposed naming:
- **Flow** = the outer unit (what we call a pipeline today)
- **Module** = the reusable composable unit (what we call a project today)
- **Step** = an instantiated module inside a flow

This keeps the meaning aligned with multi‑party routing while preserving a simple mental model for ETL‑style composition.

Notes for next spec iteration (from requirements):
- Create a new YAML spec in a `flow-spec-guide/spec/` folder with real examples modeled after existing pipeline/project YAMLs.
- Support both single-file flows (all modules inline) and split files (flow references module files).
- Allow modules to be folders (auto-discover `module.yaml` / `module.yml`) and define `module_paths`
  so flows can resolve local modules by name without hardcoding absolute paths.
- Add schema name/version (like Kubernetes) at the top for upgrade safety and forward compatibility.
- Add optional manifests (hashes) so reusable components can be resolved locally and pinned, while still allowing a looser dev mode.
- Provide native support for datasite targeting (all, subset, single), parallelism, and round-robin (e.g., `next_datasite`).
- Model start/stop conditions, waiting, and conditional completion for multi-party steps.
- Allow static inputs and runtime-dependent inputs (local environment bindings).
- Make ETL sub-steps reusable without editing upstream steps.
- Extend typing to include primitives and collections (lists, maps, sets), and file-like structures (csv, json with key paths, mappings).
- Make formats explicit so we can support singular vs batch processing (e.g., CSV row mapping, JSON key paths).
- Expose execution locations/paths (`{run}`, `{work}`, `{results}`) and support patterns like `files*.tsv` for lists.
- Protect against asset changes via hashes or manifests; allow anchors for reusable vars in YAML.
- Allow module folders with auto-discovery of `module.yaml`/`module.yml`, plus `module_paths`
  for safe local lookup by name (explicit allowlist only; no global search).

## Overlays (local patches)

We will use `.local.` overlays for on-machine overrides, similar to Kubernetes-style patching. The overlay is a separate YAML file (sidecar) that applies JSON Patch (RFC 6902) or a small strategic-merge subset. This keeps the base Flow immutable while enabling local edits.

Recommended apply order (low -> high precedence):
1. Base flow (`flow.yaml`)
2. Sidecar overlay (`flow.local.overlay.yaml`) if present
3. Runtime overlays passed in CLI order (`--overlay ...`), last wins

Example CLI:
```
bv flow run flow.yaml \
  --overlay flow.staging.overlay.yaml \
  --overlay flow.dev.overlay.yaml
```

This makes `.local.` a clear convention for local-only patches, while still allowing explicit overlays for staging/production. Runtime overlays can be inserted after the sidecar by default; if we need a different priority later, we can add a `--overlay-order` flag.

## Module discovery and safety

Recommended resolution rules for local modules:
- `source.path` may be a file or a directory. If a directory, resolve `module.yaml`/`module.yml` at its root.
- `spec.module_paths` is an explicit allowlist of local search roots. A short name like `hello` can resolve to `./modules/hello` when listed in `module_paths`.
- Resolution should be disabled unless the module ref explicitly allows local usage (e.g., `policy.allow_local: true`), to prevent accidental loading of local code.
- Avoid global recursive search; only search the configured roots and explicit paths.
