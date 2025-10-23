# Pipeline Wizard & Validator Plan

## CLI Wizard (Stage One)
- `bv pipeline create` launches an interactive flow to assemble a new `pipeline.yaml`.
- Collect pipeline metadata up front (name plus optional context literals/files).
- Declare pipeline-wide `inputs:` so required bindings (and defaults) are obvious and reusable across steps.
- Auto-populate each step's `with` block using `inputs.<name>` references and record publish mappings with type/path hints so readers see the contract immediately.
- Loop to add steps:
  - Pick a project by listing registered entries (with option to point at a filesystem path).
  - Resolve the project spec and display its declared inputs/outputs.
  - Prompt for a unique step `id` (defaulting to the project name) and allow reordering later.
  - For each input, offer available upstream outputs plus literal/glob helpers; suggest bindings when types align and warn on mismatches.
  - Let the user curate published outputs (rename, subset) while defaulting to all outputs.
  - Stub store configuration for future connectors (skippable in MVP).
- After every edit, show a validation summary so users can fix unsatisfied inputs or type conflicts before proceeding.
- Provide commands to reorder, remove, or edit steps; rerun validation after structural changes.

## Validation & Diagram Command
- Add `bv pipeline validate <pipeline.yaml>` to perform schema + type checks:
  - Ensure every referenced project exists (database lookup or relative path).
  - Confirm unique step IDs and that all required inputs are bound.
  - Apply the project type registry to flag incompatible bindings and highlight optional fallbacks.
  - Warn when nothing is published or when bindings rely on defaults only.
- Extend validation (or a dedicated `bv pipeline inspect`) to print a diagram:
  - Render an ASCII/Unicode graph of steps with project names and edge labels showing bound outputs → inputs.
  - Use styling to distinguish optional edges or type warnings.
  - Summarize published outputs and merged context sources beneath the graph for quick review.

## Implementation Notes
- Reuse existing `ProjectSpec` loading and type parsing so wizard suggestions and validator logic stay consistent with project metadata.
- Maintain an in-memory `PipelineDraft` during the wizard session, then serialize directly to `pipeline.yaml` once confirmed.
- After creation, prompt with next actions (e.g., `bv pipeline validate` followed by `bv run --set inputs.<name>=value`) to guide users toward execution once runner support lands.
- Allow steps to declare `store:` entries (starting with SQL) so the CLI can persist declared outputs into the BioVault database automatically.
