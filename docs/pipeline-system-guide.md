# BioVault Pipeline System Guide

## Overview

The BioVault pipeline system allows you to **chain multiple projects together** into multi-step workflows. Pipelines automatically:

- **Wire data flow** between steps
- **Validate type compatibility** before execution
- **Track outputs** from each step
- **Store results** in databases or files
- **Manage execution order** sequentially

A **pipeline** is defined in `pipeline.yaml` and references one or more **projects**.

---

## Pipeline Structure

### Basic Example

**pipeline.yaml:**
```yaml
name: quality-analysis-pipeline

inputs:
  samplesheet: File
  data_dir: Directory

steps:
  - id: filter
    uses: ./projects/filter-samples
    with:
      samplesheet: inputs.samplesheet
      data_dir: inputs.data_dir
    publish:
      filtered_sheet: File(filtered_results.csv)

  - id: analyze
    uses: ./projects/quality-analyzer
    with:
      samplesheet: step.filter.outputs.filtered_sheet
      data_dir: inputs.data_dir
    publish:
      analysis_results: File(results.csv)
```

**Execution:**
```bash
bv run pipeline.yaml \
  --set inputs.samplesheet=samples.csv \
  --set inputs.data_dir=./data \
  --results-dir ./results
```

---

## Pipeline Specification

### Top-Level Fields

```yaml
name: my-pipeline            # Pipeline name (required)

inputs:                      # Pipeline-level inputs (optional)
  samplesheet: File
  data_dir: Directory
  reference: File?

context:                     # Shared context (optional)
  literal:
    run_mode: production
    reviewer: analyst@lab.org
  from_json: config.json

steps:                       # Pipeline steps (required)
  - id: step1
    uses: project-reference
    with: { ... }
```

---

## Pipeline Inputs

Define reusable inputs at the pipeline level:

```yaml
inputs:
  # Simple type declaration
  samplesheet: File
  data_dir: Directory

  # Optional inputs
  config_file: File?

  # With defaults
  reference:
    type: File
    default: File(default_ref.fa)
```

**Benefits:**
- Define once, use across multiple steps
- Override at runtime with `--set inputs.<name>=<value>`
- Type validation before execution

---

## Pipeline Steps

### Step Specification

Each step in a pipeline has:

```yaml
steps:
  - id: unique_step_name       # Step identifier (required)
    uses: project-reference     # Project to run (required)
    with:                       # Input bindings (required)
      input1: value
      input2: step.other.outputs.data
    publish:                    # Output aliases (optional)
      output1: File(output.csv)
    store:                      # Database stores (optional)
      store_name:
        kind: sql
        source: output1
```

### Step ID

Unique identifier for the step:

```yaml
- id: filter_samples    # Use descriptive kebab-case names
```

Used to reference outputs: `step.filter_samples.outputs.<name>`

### Uses

Reference to a project:

```yaml
# Relative path
uses: ./projects/my-project

# Absolute path
uses: /path/to/project

# Registered project name
uses: registered-project-name
```

---

## Data Binding with `with:`

The `with:` section binds data to project inputs.

### Binding Types

#### 1. Pipeline Inputs
```yaml
with:
  samplesheet: inputs.samplesheet  # Reference pipeline input
  data_dir: inputs.data_dir
```

#### 2. Step Outputs
```yaml
with:
  samplesheet: step.filter.outputs.filtered_sheet  # Previous step output
  data_dir: inputs.data_dir
```

#### 3. Literal Values
```yaml
with:
  threshold: "0.05"                    # String literal
  min_quality: File(quality_ref.txt)   # File literal
```

#### 4. Type Constructors
```yaml
with:
  config: File(configs/default.yaml)   # File with path
  output_dir: Directory(./outputs)     # Directory path
```

### Binding Syntax

```yaml
with:
  # Format: input_name: binding_source
  samplesheet: inputs.samplesheet

  # Chain step outputs
  data: step.previous_step.outputs.result

  # Literal file path
  reference: File(/path/to/ref.fa)
```

---

## Output Publishing

Control which outputs are available downstream:

### Default Publishing
```yaml
# Without publish: all outputs available
- id: filter
  uses: ./filter-project
  with:
    data: inputs.data
  # All outputs from filter-project are available
```

### Selective Publishing
```yaml
- id: filter
  uses: ./filter-project
  with:
    data: inputs.data
  publish:
    # Only these outputs available downstream
    filtered_data: File(filtered.csv)
    metadata: File(meta.json)
```

### Renaming Outputs
```yaml
- id: analyze
  uses: ./analysis-project
  with:
    data: step.filter.outputs.filtered_data
  publish:
    # Rename output for clarity
    final_results: File(analysis_output.csv)  # Maps to 'results' output
```

**Publish Format:**
```yaml
publish:
  alias_name: Type(path)
```

- `alias_name`: Name for downstream reference
- `Type`: File, Directory, etc.
- `path`: Relative path within step results

---

## SQL Storage

Store pipeline outputs in SQLite database:

### Basic SQL Store

```yaml
- id: analyze
  uses: ./analysis-project
  with:
    data: step.filter.outputs.filtered_data
  publish:
    results: File(analysis.csv)
  store:
    results_store:              # Store name
      kind: sql                 # Store type
      destination: SQL()        # Built-in BioVault DB
      source: results           # Output to store
      table_name: analysis_{run_id}
      key_column: participant_id
```

### Store Configuration

```yaml
store:
  store_name:
    kind: sql                           # Store type (required)
    destination: SQL()                  # Database target (optional)
    source: output_name                 # Output to store (required)
    table_name: custom_table_{run_id}  # Table name (optional)
    key_column: id_column               # Primary key (optional)
    overwrite: true                     # Drop existing table (optional, default: true)
    format: csv                         # Data format (optional, auto-detected)
```

**Fields:**

- **`kind`**: Always `sql` for SQL stores
- **`destination`**: Database target
  - `SQL()` or omit → BioVault SQLite DB
  - `SQL(url:postgres://...)` → External DB (not yet supported)
- **`source`**: Name of output to store (from `publish:` or project outputs)
- **`table_name`**: Target table name
  - Supports `{run_id}` substitution
  - Defaults to `{store_name}_{run_id}`
  - Automatically prefixed with `z_results_`
- **`key_column`**: Column to use as PRIMARY KEY (optional)
- **`overwrite`**: Whether to drop existing table (default: `true`)
- **`format`**: Data format (`csv` or `tsv`, auto-detected from file extension)

### Table Naming

Tables are automatically prefixed with `z_results_`:

```yaml
table_name: analysis_{run_id}
# Creates: z_results_analysis_20251023120000
```

**Benefits:**
- Tables sort last alphabetically
- Easy to identify pipeline results
- Easy to archive/cleanup old runs

### Example with SQL Store

```yaml
name: analysis-with-storage

inputs:
  samplesheet: File
  data_dir: Directory

steps:
  - id: count_lines
    uses: ./line-counter
    with:
      samplesheet: inputs.samplesheet
      data_dir: inputs.data_dir
    publish:
      counts: File(line_counts.csv)
    store:
      line_counts_db:
        kind: sql
        destination: SQL()
        source: counts
        table_name: line_counts_{run_id}
        key_column: participant_id

  - id: summarize
    uses: ./summarizer
    with:
      counts: step.count_lines.outputs.counts
    publish:
      summary: File(summary.txt)
```

**Output:**
```
💾 Stored 'line_counts_db' output 'counts' into table z_results_line_counts_20251023120000 (rows: 150).
    source: /path/to/results/count_lines/line_counts.csv
    database: /path/to/.biovault/biovault.db
```

---

## CLI Commands

### Creating Pipelines

#### Interactive Creation
```bash
# Create pipeline with wizard
bv pipeline create

# Create with initial project
bv pipeline create --uses ./my-project --step-id first_step

# Save to custom file
bv pipeline create --file custom_pipeline.yaml
```

#### Add Steps
```bash
# Add step to existing pipeline
bv pipeline add-step

# Add to specific file
bv pipeline add-step --file my_pipeline.yaml
```

### Validating Pipelines

```bash
# Basic validation
bv pipeline validate pipeline.yaml

# Validation with diagram
bv pipeline validate --diagram pipeline.yaml
```

**Validation Checks:**
- ✅ All referenced projects exist
- ✅ Input bindings match types
- ✅ Required inputs are satisfied
- ✅ Step IDs are unique
- ✅ Output references are valid
- ✅ SQL stores reference valid outputs

### Running Pipelines

```bash
# Basic run
bv run pipeline.yaml \
  --set inputs.samplesheet=data.csv \
  --set inputs.data_dir=./data

# Custom results directory
bv run pipeline.yaml \
  --set inputs.samplesheet=data.csv \
  --results-dir ./custom_results

# Override step inputs
bv run pipeline.yaml \
  --set inputs.samplesheet=data.csv \
  --set filter.threshold=0.01

# Dry run
bv run pipeline.yaml --dry-run \
  --set inputs.samplesheet=data.csv
```

### Inspecting Pipelines

```bash
# View pipeline structure
bv pipeline inspect pipeline.yaml

# Detailed validation output
bv pipeline validate --diagram pipeline.yaml
```

---

## Execution Model

### Sequential Execution

Steps run in order, one at a time:

```yaml
steps:
  - id: step1      # Runs first
    uses: project-a

  - id: step2      # Runs after step1 completes
    uses: project-b
    with:
      data: step.step1.outputs.result

  - id: step3      # Runs after step2 completes
    uses: project-c
```

### Results Directory Structure

```
results/
├── pipeline-name/
│   ├── step1/
│   │   ├── output1.csv
│   │   └── output2.txt
│   ├── step2/
│   │   └── result.csv
│   └── step3/
│       └── final.txt
```

Or with custom `--results-dir`:

```
custom_results/
├── step1/
├── step2/
└── step3/
```

---

## Advanced Features

### Pipeline Context

Share configuration across all steps:

```yaml
context:
  literal:
    run_mode: production
    quality_threshold: 30
    reviewer: analyst@example.com

  from_json: shared_config.json  # Load from JSON file
```

**Access in projects:**
```groovy
workflow USER {
    take:
        context
        // ...

    main:
        def runMode = context.run_mode
        def threshold = context.quality_threshold
}
```

### Optional Inputs

Mark inputs as optional:

```yaml
inputs:
  required_data: File
  optional_config: File?      # Optional
  optional_dir: Directory?    # Optional
```

### Type Validation

Pipeline validates compatibility:

```yaml
# ✅ Valid: File → File
- id: step1
  outputs:
    data: File

- id: step2
  with:
    input: step.step1.outputs.data  # File → File ✅

# ❌ Invalid: File → Directory
- id: step3
  with:
    input_dir: step.step1.outputs.data  # File → Directory ❌
```

---

## Complete Example

**pipeline.yaml:**
```yaml
name: variant-analysis-pipeline

inputs:
  samplesheet: File
  data_dir: Directory
  reference_genome: File?

steps:
  - id: filter_quality
    uses: ./projects/quality-filter
    with:
      samplesheet: inputs.samplesheet
      data_dir: inputs.data_dir
    publish:
      filtered_samples: File(filtered.csv)

  - id: count_variants
    uses: ./projects/variant-counter
    with:
      samplesheet: step.filter_quality.outputs.filtered_samples
      data_dir: inputs.data_dir
      reference: inputs.reference_genome
    publish:
      variant_counts: File(counts.csv)
    store:
      variant_counts_db:
        kind: sql
        destination: SQL()
        source: variant_counts
        table_name: variants_{run_id}
        key_column: sample_id

  - id: summarize
    uses: ./projects/summary-stats
    with:
      counts: step.count_variants.outputs.variant_counts
    publish:
      summary_report: File(summary.txt)
      stats_table: File(stats.csv)
    store:
      summary_stats_db:
        kind: sql
        source: stats_table
        table_name: summary_stats_{run_id}
```

**Running:**
```bash
# Run pipeline with required inputs
bv run pipeline.yaml \
  --set inputs.samplesheet=samples.csv \
  --set inputs.data_dir=./sequencing_data \
  --results-dir ./analysis_results

# Run with optional reference
bv run pipeline.yaml \
  --set inputs.samplesheet=samples.csv \
  --set inputs.data_dir=./data \
  --set inputs.reference_genome=hg38.fa \
  --results-dir ./results
```

**Validation:**
```bash
bv pipeline validate --diagram pipeline.yaml
```

**Output:**
```
🔑 Pipeline inputs:
  - samplesheet : File
  - data_dir : Directory
  - reference_genome : File?

📦 Steps:
  • filter_quality → quality-filter (./projects/quality-filter)
      samplesheet ← inputs.samplesheet [ok] File
      data_dir ← inputs.data_dir [ok] Directory

  • count_variants → variant-counter (./projects/variant-counter)
      samplesheet ← step.filter_quality.outputs.filtered_samples [ok] File → File
      data_dir ← inputs.data_dir [ok] Directory
      reference ← inputs.reference_genome [ok] File?
      store:
        variant_counts_db → sql(table: variants_{run_id}, source: variant_counts)

  • summarize → summary-stats (./projects/summary-stats)
      counts ← step.count_variants.outputs.variant_counts [ok] File → File
      store:
        summary_stats_db → sql(table: summary_stats_{run_id}, source: stats_table)

✅ Pipeline is valid
```

---

## Best Practices

### 1. Descriptive Step IDs
```yaml
# Good
- id: filter_low_quality_samples
- id: count_variants_per_sample
- id: generate_summary_statistics

# Avoid
- id: step1
- id: process
- id: analyze
```

### 2. Explicit Input Bindings
```yaml
# Good - clear data flow
- id: analyze
  with:
    input_data: step.filter.outputs.filtered_samples
    reference: inputs.reference_genome

# Avoid - unclear source
- id: analyze
  with:
    data: step.filter.outputs.result  # Which result?
```

### 3. Meaningful Output Aliases
```yaml
# Good
publish:
  high_quality_variants: File(filtered_variants.vcf)
  variant_summary_stats: File(summary.csv)

# Avoid
publish:
  output1: File(out1.txt)
  result: File(data.csv)
```

### 4. Use Pipeline Inputs for Reusability
```yaml
# Good - reusable across steps
inputs:
  reference_genome: File
  quality_threshold: String

steps:
  - id: step1
    with:
      reference: inputs.reference_genome
  - id: step2
    with:
      reference: inputs.reference_genome  # Reuse

# Avoid - duplicated literal paths
steps:
  - id: step1
    with:
      reference: File(/path/to/ref.fa)
  - id: step2
    with:
      reference: File(/path/to/ref.fa)  # Duplicated
```

### 5. SQL Store Naming
```yaml
# Good - descriptive, uses run_id
store:
  variant_analysis_results:
    table_name: variant_counts_{run_id}
    key_column: sample_id

# Avoid - static names (overwrites previous runs)
store:
  results:
    table_name: results  # No run_id!
```

---

## Troubleshooting

### Common Issues

**1. Type Mismatch**
```
Error: Type mismatch: expected File but found Directory
```
→ Check binding types match between steps

**2. Missing Input**
```
Error: Pipeline input 'samplesheet' is not set
```
→ Provide value with `--set inputs.samplesheet=<path>`

**3. Unknown Output Reference**
```
Error: Step 'step1' output 'data' not found
```
→ Check `publish:` section of referenced step

**4. Step Without Project**
```
Error: Step 'analyze' has no project to run
```
→ Add `uses:` field to step or register the project

**5. SQL Store Source Not Found**
```
Error: Store 'db_store' references unknown output 'results'
```
→ Ensure `source:` matches name in `publish:` or project outputs

---

## Next Steps

- Learn about [Project System](project-system-guide.md) to create reusable components
- See [YAML Types Reference](yaml-types-reference.md) for type specifications
- Check [Pipeline Design Patterns](nextflow_project_pipeline_design.md) for advanced usage
