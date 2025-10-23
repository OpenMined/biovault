# BioVault Project System Guide

## Overview

The BioVault project system provides a structured way to define reusable computational workflows using Nextflow. A **project** is a self-contained unit that declares:

- **Inputs** it expects (files, directories, parameters)
- **Outputs** it produces
- **Assets** it uses (scripts, schemas, reference data)
- **Workflow logic** written in Nextflow DSL2

Projects are portable, testable, and composable—you can chain them together into pipelines.

---

## Project Structure

A typical project directory contains:

```
my-project/
├── project.yaml      # Project specification (required)
├── workflow.nf       # Nextflow workflow logic (required)
├── assets/           # Scripts, schemas, reference files
│   ├── process.py
│   └── schema.yaml
└── README.md         # Documentation (optional)
```

### Core Files

#### 1. `project.yaml`

The project specification defines the contract for your workflow:

```yaml
name: my-analysis-project
author: user@example.com
workflow: workflow.nf
template: dynamic-nextflow
version: 1.0.0

assets:
  - process.py
  - schema.yaml

parameters:
  - name: threshold
    type: String
    description: Filtering threshold value
    default: "0.05"

inputs:
  - name: samplesheet
    type: File
    description: CSV file with sample metadata
    format: csv

  - name: data_dir
    type: Directory
    description: Directory containing sample data files

outputs:
  - name: results
    type: File
    description: Analysis results CSV
    format: csv
    path: analysis_results.csv

  - name: summary
    type: File
    description: Summary statistics
    path: summary.txt
```

**Key Sections:**

- **`name`**: Unique identifier for this project
- **`author`**: Contact email for the project maintainer
- **`workflow`**: Path to the Nextflow workflow file (usually `workflow.nf`)
- **`template`**: Set to `dynamic-nextflow` to use BioVault's template system
- **`version`**: Semantic version (e.g., `1.0.0`)
- **`assets`**: List of files to include (scripts, configs, reference data)
- **`parameters`**: User-configurable options with defaults
- **`inputs`**: Data inputs required by the workflow
- **`outputs`**: Data outputs produced by the workflow

#### 2. `workflow.nf`

The Nextflow workflow implements your computational logic:

```groovy
nextflow.enable.dsl=2

workflow USER {
    take:
        context       // BiovaultContext - auto-injected
        samplesheet   // File - matches inputs section
        data_dir      // Directory - matches inputs section

    main:
        // Access parameters via context
        def assetsDir = file(context.params.assets_dir)
        def threshold = context.params.threshold ?: '0.05'

        // Your workflow logic here
        results_ch = process_samples(samplesheet, data_dir, threshold)

    emit:
        results = results_ch    // Must match outputs in project.yaml
        summary = summary_ch
}

process process_samples {
    publishDir params.results_dir, mode: 'copy'

    input:
        file sheet
        val data_dir
        val threshold

    output:
        path 'analysis_results.csv'

    script:
    """
    python3 ${assetsDir}/process.py \
        --input ${sheet} \
        --data-dir ${data_dir} \
        --threshold ${threshold} \
        --output analysis_results.csv
    """
}
```

**Important Conventions:**

1. **Workflow must be named `USER`** - BioVault generates wrapper code that calls this
2. **First parameter is always `context`** - Auto-injected BioVaultContext
3. **Remaining parameters match `inputs` in order** - Defined in `project.yaml`
4. **Emit outputs by name** - Names must match `outputs` section

---

## CLI Commands

### Creating Projects

#### Interactive Creation
```bash
# Create a new project interactively
bv project create

# Create from a specific directory
bv project create --path ./my-project
```

#### List Examples
```bash
# See available project templates
bv project examples
```

### Managing Projects

#### Register a Project
```bash
# Register local project for reuse
bv project import ./path/to/project

# Register with custom name
bv project import ./path/to/project --name my-custom-name
```

#### List Projects
```bash
# List all registered projects
bv project list
```

#### View Project Details
```bash
# View project specification
bv project view ./my-project

# Show detailed project info
bv project show my-project-name
```

#### Delete Project
```bash
# Remove a registered project
bv project delete my-project-name
```

### Running Projects

#### Basic Run
```bash
# Run a project with inputs
bv run ./my-project \
  --samplesheet data/samples.csv \
  --data_dir data/files

# Run with custom parameters
bv run ./my-project \
  --samplesheet data/samples.csv \
  --data_dir data/files \
  --param.threshold 0.01
```

#### Advanced Options
```bash
# Custom results directory
bv run ./my-project \
  --samplesheet data.csv \
  --results-dir custom_results

# Dry run (show commands without executing)
bv run ./my-project --dry-run --samplesheet data.csv

# Resume from previous run
bv run ./my-project --resume --samplesheet data.csv
```

---

## Input & Output Types

### Supported Types

| Type | Description | Example |
|------|-------------|---------|
| `File` | Single file | `data.csv`, `config.json` |
| `Directory` | Directory path | `/path/to/data/` |
| `String` | Text parameter | `"threshold_value"` |
| `Bool` | Boolean flag | `true`, `false` |

### Optional Types

Add `?` to make any type optional:

```yaml
inputs:
  - name: optional_config
    type: File?
    description: Optional configuration file
```

### Type Formats

For `File` types, specify format for better validation:

```yaml
outputs:
  - name: results
    type: File
    format: csv     # csv, tsv, json, txt, etc.
    path: output.csv
```

---

## Parameters vs Inputs

### Parameters
- **User-configurable settings** (thresholds, flags, options)
- Have **default values**
- Accessed via `context.params.<name>`
- Examples: `threshold`, `enable_qc`, `output_format`

```yaml
parameters:
  - name: min_quality
    type: String
    default: "30"
    description: Minimum quality score
```

### Inputs
- **Required data** (files, directories)
- **No defaults** (must be provided at runtime)
- Passed as workflow parameters
- Examples: `samplesheet`, `data_dir`, `reference_genome`

```yaml
inputs:
  - name: samplesheet
    type: File
    description: Sample metadata CSV
    format: csv
```

---

## Assets

Assets are files bundled with your project:

```yaml
assets:
  - process_data.py      # Python script
  - filter.R             # R script
  - schema.yaml          # Data schema
  - reference.fa         # Reference file
```

**Accessing Assets in Workflow:**

```groovy
def assetsDir = file(context.params.assets_dir)
def script = file("${assetsDir}/process_data.py")
```

Assets are automatically copied to the execution environment.

---

## BiovaultContext

Every workflow receives a `BiovaultContext` as the first parameter:

```groovy
workflow USER {
    take:
        context  // BiovaultContext - always first
        // ... other inputs

    main:
        // Access context fields
        def runId = context.run_id
        def datasite = context.datasite
        def assetsDir = file(context.params.assets_dir)
        def resultsDir = file(params.results_dir)

        // Access user parameters
        def threshold = context.params.threshold
}
```

**Context Fields:**

- `run_id`: Unique run identifier
- `datasite`: Current datasite name
- `user`: User email
- `run_timestamp`: Run start time
- `params`: Map of all parameters
- `assets_dir`: Path to assets directory

---

## Best Practices

### 1. Clear Naming
```yaml
name: descriptive-project-name  # Use kebab-case
author: your-email@example.com
```

### 2. Comprehensive Descriptions
```yaml
inputs:
  - name: samplesheet
    type: File
    description: CSV with columns: participant_id, file_path, metadata
    format: csv
```

### 3. Meaningful Parameters
```yaml
parameters:
  - name: quality_threshold
    type: String
    default: "30"
    description: Minimum phred quality score (0-40)
```

### 4. Explicit Output Paths
```yaml
outputs:
  - name: filtered_data
    type: File
    path: filtered_results.csv  # Explicit filename
    format: csv
```

### 5. Version Your Projects
```yaml
version: 1.2.0  # major.minor.patch
```

---

## Example: Complete Project

**project.yaml:**
```yaml
name: sample-quality-filter
author: analyst@lab.org
workflow: workflow.nf
template: dynamic-nextflow
version: 1.0.0

assets:
  - filter_quality.py

parameters:
  - name: min_quality
    type: String
    default: "30"
    description: Minimum quality score threshold

inputs:
  - name: input_data
    type: File
    description: Raw sequencing data
    format: csv

outputs:
  - name: filtered_data
    type: File
    description: Quality-filtered results
    format: csv
    path: filtered.csv
```

**workflow.nf:**
```groovy
nextflow.enable.dsl=2

workflow USER {
    take:
        context
        input_data

    main:
        def assetsDir = file(context.params.assets_dir)
        def filterScript = file("${assetsDir}/filter_quality.py")
        def minQuality = context.params.min_quality

        filtered_ch = filter_quality(
            input_data,
            filterScript,
            Channel.value(minQuality)
        )

    emit:
        filtered_data = filtered_ch
}

process filter_quality {
    publishDir params.results_dir, mode: 'copy'

    input:
        path input_file
        path filter_script
        val min_qual

    output:
        path 'filtered.csv'

    script:
    """
    python3 ${filter_script} \
        --input ${input_file} \
        --min-quality ${min_qual} \
        --output filtered.csv
    """
}
```

**Running:**
```bash
bv run ./sample-quality-filter \
  --input_data raw_data.csv \
  --param.min_quality 35
```

---

## Troubleshooting

### Common Issues

**1. Workflow outputs don't match project.yaml**
```
Error: Output 'results' not emitted by workflow
```
→ Ensure `emit:` block matches `outputs:` names exactly

**2. Assets not found**
```
Error: File not found: process.py
```
→ List all files in `assets:` section of project.yaml

**3. Parameter not accessible**
```
Error: No such property: threshold
```
→ Access via `context.params.threshold`, not just `threshold`

**4. Type mismatch**
```
Error: Expected File, got Directory
```
→ Check input types in project.yaml match CLI arguments

---

## Next Steps

- Learn about [Pipeline System](pipeline-system-guide.md) to chain projects
- See [YAML Types Reference](yaml-types-reference.md) for detailed type specifications
- Check [Nextflow Project Design](nextflow_project_pipeline_design.md) for advanced patterns
