# Nextflow Project & Pipeline DSL Design

## Goals
- Allow each project to describe its runtime contract (inputs, outputs, assets) without hard-coded templates.
- Autogenerate the outer Nextflow scaffolding that maps declared types onto channels/values so inner workflows stay minimal.
- Enable composition of multiple projects into multi-step pipelines with static validation of type compatibility.
- Keep authoring experience friendly: interactive wizard plus JSON/YAML driven creation for automation.

## Current Pain Points
- Templates in `cli/src/templates` are static; adding new variants requires editing Groovy by hand.
- Inputs/outputs are implicit (e.g. sheet template always expects a CSV) which makes chaining projects brittle.
- No shared language for domain types (ParticipantSheet, GenotypeRecord, etc.), so downstream steps cannot reason about upstream outputs.
- Linking projects today means manual channel wiring inside ad-hoc Nextflow code.

## Project Specification (`project.yaml`)
Top-level keys stay familiar, with `template` switching to a dynamic generator.

```yaml
name: apol1-sheet
author: madhava@openmined.org
workflow: workflow.nf
template: dynamic-nextflow
version: 1.0.0
assets:
  - schema.yaml
  - apol1_classifier.py

parameters:
  - name: hereditary_panel
    type: Bool
    default: true
    description: Run hereditary breast cancer panel

inputs:
  - name: rows
    type: List[GenotypeRecord]
    description: List of genotype records to score

outputs:
  - name: scored_sheet
    type: ParticipantSheet
    format: csv
    path: results/scored.csv
    description: Participant sheet including classifier output

  - name: scored_assets
    type: Directory
    path: results/artifacts
    description: Generated artifacts for downstream reuse
```

The optional `version` field follows semver and lets us track template expectations for upgrades/back-compat.

Parameter entries declare user-facing toggles or scalars. The generator adds them to the Biovault context and exposes them as plain values (strings/bools); enum selections resolve to strings.

`BiovaultContext` is injected automatically into the generated wrapper; project authors do not list it under `inputs` unless they want to rename it in their workflow signature.

### Parameter fields
- `name` (required): identifier shown to users and surfaced in the context as `params.<name>`.
- `type` (required): currently `String`, `Bool`, or `Enum[...]`.
- `description` (optional): guidance for UI/wizard prompts.
- `default` (optional): value used when not supplied at runtime.
- `choices` (optional, only for enums): allowed literals.
- `advanced` (optional Bool): hide from simple wizards unless the user opts in.

### Input object fields
- `name` (required): identifier exposed to workflow `take:` block.
- `type` (required): Python-expressible type signature (see Type System below).
- `description` (optional): documentation for wizard/help output.
- Inputs are bound from the CLI or pipeline by name; the `BiovaultContext` is auto-injected without being listed.
- To expose paths like `assets_dir` or `results_dir` to your workflow, add them as explicit inputs (e.g. `Directory` or `String` types); they are no longer bundled into the context automatically.
- `format` (optional): hint for loaders (`csv`, `tsv`, `json`, `directory`, `binary`).
- `path` / `pattern` (optional): used when loaders expect a filesystem binding (e.g. directory glob).
- `default` (optional): literal default when absent (honours type coercion).

### Output object fields
- `name` (required): key the pipeline can reference (e.g. `step.outputs.scored_sheet`).
- `type` (required): same type language as inputs.
- `description` (optional).
- `format` (optional) and `path` (optional) describe where artifacts land relative to `results_dir`.
- `visibility` (optional): `public`, `private`, `internal` to control sharing and hashing.

Outputs default to being written under `${results_dir}/${name}` when `path` is not supplied; the generator will create directories as needed.

## Type System
We treat types as declarative descriptors resolved by a registry in the CLI. Categories:

Note: append `?` to any type to declare it optional (e.g., `String?`, `List[File]?`).

| Type | Nextflow binding | Notes |
|------|------------------|-------|
| `String` | Plain value passed into `USER` | For parameters and literal flags (includes enum selections). |
| `Bool` | Plain value passed into `USER` | For feature toggles surfaced via `parameters`. |
| `File` | `Channel.fromPath` emitting `path` objects | Read-only staging; use `publishDir`/results for writes. |
| `Directory` | `Channel.fromPath` emitting directory `path` objects | Exposed as channels so steps can fan-out listings. |
| `ParticipantSheet` | Channel emitting row maps plus sheet metadata | Built via sheet loader; also exposes lazy index by `participant_id`. |
| `GenotypeRecord` | Channel emitting immutable maps with required keys | `genotype_file` is a staged `path`; helper index keyed by `participant_id`. |
| `BiovaultContext` | `Channel.value` of deep-frozen map | Runtime metadata only; directories provided separately if declared. |
| `List[T]` | Queue channel streaming each `T` | Use for collections; generator expands loaders accordingly. |
| `Map[String, T]` | `Channel.value` of immutable map | Handy for lookup tables (e.g., participant → record). |
| `T?` | Follows binding of `T` but allows null/empty | Defaults respected; channels may emit zero items. |

All non-`String` bindings flow through channels so processes can stream data. Composite records (e.g., `GenotypeRecord`) keep their file handles as `path` objects inside the immutable map, preserving channel semantics while still offering dictionary-style lookups.

### Domain structures
- `GenotypeRecord` (record type)
  - Represented as immutable map keyed by the column names below.
  - `participant_id: String`
  - `genotype_file: File` (Nextflow `path` item)
  - `grch_build: String?`
  - `source: String?`

- `ParticipantSheet` (dataset)
  - Physical representation: CSV/TSV with schema metadata.
  - Logical representation: immutable map with keys `rows` (channel of row maps), `header`, and `index` (lookup by `participant_id`).
  - Convenience conversion `ParticipantSheet.rows -> List[GenotypeRecord]` when genotype columns are available; generator raises a validation warning otherwise.

- `BiovaultContext`
  - Always injected even if not declared; explicit declaration gives typed access in workflow.
  - Keys: `run_id`, `datasite`, `user`, `run_timestamp`, auth hints, user-provided overrides from the pipeline/root CLI, and `params` (map of declared project parameters), plus other deployment metadata (no asset or results paths).
  - Delivered as a deep-frozen map (`asImmutable()` applied recursively) so inner workflows can read but not mutate the shared context.

### Type registry responsibilities
The CLI maintains a registry describing:
- How to parse CLI/pipeline inputs into Nextflow params.
- Snippets to emit in `template.nf` for each type (`fromPath`, CSV parsing, JSON decoding).
- Compatibility graph (e.g. `ParticipantSheet` is assignable to `List[GenotypeRecord]`).
- Parameters -> context hydration (apply defaults, CLI overrides, literal/JSON merges).
- Validation hooks executed during `bv project verify` and `bv pipeline validate`.

## Dynamic Template Generation
When `template: dynamic-nextflow`, `bv project create` (or `bv project build`) produces `template.nf` based on the spec.

Generation steps:
1. Read `project.yaml`, resolve types using registry, and expand imports required by loaders.
2. Declare `params.<name>` only for the inputs/outputs actually referenced (plus internal plumbing for outputs).
3. Emit loader code per input:
   - Scalars -> `Channel.value(params.<name>)` with optional casting.
   - Files/Directories -> `Channel.fromPath` with existence checks.
   - Sheets -> reuse the existing CSV loader (currently in sheet template) moved into a reusable library module.
   - Context -> build map using runtime metadata (datasite, run id, etc.) without implicitly adding asset/results paths.
4. Emit any helper includes (e.g. `include { SheetLoader } from "${lib}").
5. Construct `workflow USER { take: ... }` block so the invocation list is `[context] + declared inputs` in order.
6. Route outputs: create placeholders or helper channels so the user workflow can `emit` using the declared output identifiers.
7. Auto-create `results_dir` internally for file staging, but do not pass it unless the project declares an input for it.

### Example generated scaffold (excerpt)
```nextflow
nextflow.enable.dsl=2
include { SheetLoader } from params._bv_lib

params.results_dir = params.results_dir ?: 'results'
params.rows = params.rows ?: null
params.context_json = params.context_json ?: null

results_dir = file(params.results_dir)
results_dir.mkdirs()

rows_ch = SheetLoader.load(params.rows, schema: file('${params.assets_dir}/schema.yaml'))
context_ch = Channel.value( SheetLoader.buildContext(params.context_json).asImmutable() )

include { USER } from "${params.work_flow_file}"

workflow {
    USER(
        context_ch,
        rows_ch
    )
}
```
The generated file contains only wiring; project authors implement `workflow.nf` without repeating boilerplate.

## Outputs in Nextflow
Outputs are declared so downstream tooling knows what to collect. The generator provides helpers:
- Injects a Groovy helper `emit_<name>(path)` or expects the workflow to emit channels with names matching `outputs`.
- Adds validation step at the end of the workflow to ensure declared paths exist (unless marked optional).
- Annotates metadata (type, format) in a `.bv/project_state.json` for the CLI to read when linking steps.
- Optionally runs a post-check that no files outside `results_dir` were touched, letting us flag accidental writes back into assets or working directories.

## Pipeline Composition (`pipeline.yaml`)
A separate file describes orchestration of multiple projects.

Each referenced project input is bound by name. The shared `BiovaultContext` arrives automatically, so the pipeline only wires explicit data dependencies. A step can optionally curate its outputs via `publish`, keeping the rest internal.

```yaml
name: apol1-with-stats
workdir: work/pipelines/apol1

inputs:
  rows: File
  data_root: Directory

context:
  literal:
    run_mode: production
    reviewer: madhava@openmined.org
  from_json: configs/shared_context.json

steps:
  - id: preprocess
    uses: preprocess
    where: datasite-default
    with:
      files: FileGlob("/data/input/*.vcf")
    publish:
      sheet: filtered_sheet
      files: filtered_files

  - id: sheet
    uses: apol1-sheet
    with:
      rows: step.preprocess.outputs.sheet
      data_dir: inputs.data_root
    publish:
      scored: scored_sheet

  - id: stats
    uses: stats
    with:
      rows: step.sheet.outputs.scored
    publish:
      stats: stats_txt
    store:
      stats_sql:
        kind: sql
        destination: SQL()
        source: stats
        table_name: stats_{run_id}
        key_column: participant_id

complete:
  expect:
    - step.stats.outputs.stats
```

`context.literal` holds an inline YAML map of extra arguments merged into the Biovault context, while `context.from_json` (optional) points at a file that will be parsed and merged on top (JSON or YAML). CLI-provided values win last, so runs can override pipeline defaults without editing the spec.

`store` entries describe side-effects to perform once a step succeeds. In the example above, the `stats` output is written to a local SQLite database managed by Biovault, with `key_column` pointing to the participant column used for upserts. Stores can target other connectors (object storage, APIs) in future revisions.

The CLI now supports the `sql` store kind for local BioVault data. During `bv run pipeline.yaml`, CSV/TSV outputs referenced by a `sql` store are loaded and written into the BioVault SQLite database (default table name supports `{run_id}` substitution). Pipelines can pass overrides such as `--set inputs.samplesheet=/path` and the runner handles ingestion automatically.

### Pipeline schema highlights
- `context`: optional block providing defaults for the shared `BiovaultContext` via `literal` maps or `from_json` file references.
- `inputs`: optional block declaring reusable bindings (`File`, `Directory`, etc.) that steps can reference as `inputs.<name>`.
- `steps[].id`: unique handle for wiring.
- `steps[].uses`: reference to a project directory or registry entry.
- `steps[].with`: map of input bindings. Values can be literals, references (`step.<id>.outputs.<name>`), or helper constructors (`FileGlob`, `StaticPipe`, `Literal`).
- `steps[].publish`: optional map selecting which outputs to expose and how to rename them for downstream steps (defaults to all outputs under their declared names).
- `steps[].store`: optional map describing side-effect destinations (e.g., SQL export) executed after a step completes.
- `steps[].where`: optional execution target (datasite, queue, tenant) for future multi-site orchestration.
- Control flow additions for later phases: `when` (conditional), `foreach` (iterate over list), `scatter`/`gather` semantics.
- `pipeline.context`: shared `BiovaultContext` seeded from CLI flags (datasite, submission id, etc.), optionally augmented with literal overrides or JSON files, and auto-injected into every step (no wiring needed).

When `publish` is omitted, every declared project output is available as `step.<id>.outputs.<name>`. Providing a `publish` map lets authors curate/rename the subset they want downstream (the files still exist on disk even if not published).
To keep that mapping easy to review without encouraging edits, the CLI can regenerate a read-only summary (e.g. `pipeline.lock.yaml` or `bv pipeline inspect`) showing each `step.<id>.outputs.<key>` binding after publish rules are applied. The pipeline spec remains the single source of truth.

Validation ensures:
- Every referenced project exists and its spec is parseable.
- Each binding's type is compatible (using registry conversions).
- Required inputs are satisfied; optional ones may fall back to defaults.
- For `List` outputs feeding `List` inputs, cardinality remains `many`. For mismatches (e.g. `ParticipantSheet` into `List[GenotypeRecord]`) the registry inserts a conversion node if possible, otherwise surfaces an error with remediation hints.

## Biovault SQL Integration (Planned)
- Introduce an `SQL` loader that can hydrate `ParticipantSheet` or `List[GenotypeRecord]` inputs directly from the Biovault database (e.g., `SQL(location: biovault.sqlite, query: "select * from genotype_view")`).
- Query results stream through channels the same way filesystem-backed inputs do; large result sets remain back-pressure-friendly.
- Companion `store` handlers persist outputs back to SQL. `key_column` identifies which column becomes the primary key (defaults to `participant_id` when present); additional columns come from the emitted sheet/records.
- Runs write into run-scoped tables (e.g., `stats_{run_id}`) while also supporting merging into canonical tables via configuration.
- The CLI will manage schema creation/migration and capture metadata so downstream steps (or external consumers) can look up exactly which run produced a dataset.

## CLI Experience
### Project wizard (`bv project create`)
1. Collect name, author (prefill from config), workflow filename.
2. Prompt for parameters (name, type, default, description) to surface runtime toggles.
3. Prompt for inputs one by one (name, choose type from registry shortlist, description, defaults).
4. Prompt for outputs similarly (name, type, destination hints).
5. Confirm assets to copy.
6. Emit `project.yaml`, `workflow.nf` stub with `take:` signature, `template.nf` generated from spec.
7. Optional `--spec project.json` to bypass prompts; JSON matches YAML schema.

### Supporting commands
- `bv project inspect <dir>`: show inputs/outputs/types rendered nicely.
- `bv project build <dir>`: regenerate `template.nf` after manual edits to spec.
- `bv pipeline validate <pipeline.yaml>`: type-check bindings, ensure referenced outputs exist.
- `bv run <pipeline.yaml> --set inputs.<name>=<value>`: run pipelines sequentially, wiring shared inputs and optional `store` targets (e.g., SQL exports) automatically.

## Implementation Phases
1. **Registry groundwork**: define Rust structs for types, loaders, compatibility rules; port sheet loader helpers into reusable module.
2. **Spec parser & validator**: load extended `project.yaml`, enforce schema, emit helpful errors.
3. **Template generator**: produce `template.nf` + updated `workflow.nf` stub using registry metadata.
4. **Output metadata capture**: standardise where each run records produced artifacts (likely `.bv/output.json`).
5. **Pipeline DSL MVP**: parse `pipeline.yaml`, honour default publish behaviour, execute steps sequentially on a single runtime.
6. **CLI UX**: wizard prompts, `--spec` support, pretty printers.
7. **Store connectors**: implement SQL (and future) backends for `store`, including schema management and run metadata tracking.
8. **Advanced orchestration** (future): remote datasite execution, parallel scatters, completion conditions.

## Open Questions
- How opinionated should loaders be? e.g. `List[File]` from directory glob vs explicit array of paths supplied by user.
- Should we auto-materialise `ParticipantSheet -> List[GenotypeRecord]`, or require explicit conversion step to keep data transformations visible?
- Where do shared Groovy helpers (`SheetLoader`, `TypeCoercions`) live: vendored per project or referenced from the CLI install?
- How do we version the type registry so existing projects remain reproducible? (Proposal: stamp generator version in `project.yaml` and keep backward-compatible loaders.)
- What is the exact shape of the SQL connector configuration (table naming, column typing, transaction safety) so it plays nicely with Biovault migrations?
- Do we provide a generated participant→record index automatically, or require workflows to request the `Map[String, GenotypeRecord]` explicitly to signal lookups?
- Output discovery: do we enforce declared `path` to match actual `emit`, or allow runtime to override (with metadata update)?
- What is the minimal context contract we guarantee across deployments (local, datasite, remote environments)?

## Next Steps for Review
- Align on initial type list (especially domain-specific ones) and compatibility rules.
- Decide default behaviour for `List[File]` and `Directory` inputs (single Channel vs map of staged files).
- Prototype generator on the three sample projects (`preprocess`, `apol1-sheet`, `stats`) to smoke test ergonomics.
- Iterate on pipeline syntax once individual project scaffolding feels solid.
