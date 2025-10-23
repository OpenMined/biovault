# BioVault YAML Types & Syntax Reference

## Overview

BioVault uses a type system to bridge YAML configuration files, filesystem paths, and Nextflow Groovy code. This guide explains:

- YAML type syntax for projects and pipelines
- How types map to files and directories
- How types convert to Nextflow channels and values
- Type compatibility rules

---

## Type System Basics

### Core Principles

1. **Types describe data shape** - File, Directory, String, etc.
2. **Types enforce compatibility** - Pipeline validation checks type matches
3. **Types map to Nextflow** - Automatic conversion to channels
4. **Types are declarative** - No runtime type casting needed

### Type Categories

| Category | Types | Use Case |
|----------|-------|----------|
| **Scalar** | `String`, `Bool` | Configuration, parameters |
| **Filesystem** | `File`, `Directory` | Data files, folders |
| **Optional** | `Type?` | Optional inputs/outputs |

---

## Type Definitions

### 1. File

**Description:** Single file on filesystem

**YAML Syntax:**
```yaml
# In project.yaml inputs/outputs
- name: samplesheet
  type: File
  format: csv          # Optional format hint
  path: output.csv     # Optional path (outputs only)
```

**Pipeline Binding:**
```yaml
# Literal file
input_file: File(data/samples.csv)

# Reference pipeline input
input_file: inputs.samplesheet

# Reference step output
input_file: step.filter.outputs.filtered_data
```

**Filesystem Mapping:**
- Points to **single file**
- Absolute or relative path
- Must exist (for inputs) or be created (for outputs)

**Nextflow Conversion:**
```groovy
// In workflow.nf - File becomes path object
workflow USER {
    take:
        context
        input_file   // Nextflow path object

    main:
        // Use directly in processes
        process_file(input_file)
}

process process_file {
    input:
        path input_csv    // Nextflow stages the file

    output:
        path 'result.csv'

    script:
    """
    process.py ${input_csv}
    """
}
```

**CLI Usage:**
```bash
# Pass file path
bv run ./project --input_file data.csv

# Pipeline override
bv run pipeline.yaml --set inputs.samplesheet=samples.csv
```

---

### 2. Directory

**Description:** Directory/folder on filesystem

**YAML Syntax:**
```yaml
- name: data_dir
  type: Directory
  description: Folder containing sample files
```

**Pipeline Binding:**
```yaml
# Literal directory
data_dir: Directory(/path/to/data)

# Reference
data_dir: inputs.data_directory
```

**Filesystem Mapping:**
- Points to **directory path**
- Contains multiple files
- Directory structure preserved

**Nextflow Conversion:**
```groovy
workflow USER {
    take:
        context
        data_dir     // String path to directory

    main:
        // Pass as value (directory path)
        def dataDirPath = data_dir
        process_directory(Channel.value(dataDirPath))

        // Or use for file discovery
        def files = file("${data_dir}/*.txt")
}

process process_directory {
    input:
        val dir_path

    script:
    """
    ls ${dir_path}/*.csv > file_list.txt
    """
}
```

**CLI Usage:**
```bash
bv run ./project --data_dir ./sequencing_data
```

---

### 3. String

**Description:** Text value (parameters, configuration)

**YAML Syntax:**
```yaml
parameters:
  - name: threshold
    type: String
    default: "0.05"
    description: Quality threshold value
```

**Pipeline Binding:**
```yaml
# Literal string
threshold: "0.05"

# Reference parameter (via context)
# Accessed in workflow as: context.params.threshold
```

**Nextflow Conversion:**
```groovy
workflow USER {
    take:
        context

    main:
        // Access string parameters via context
        def threshold = context.params.threshold ?: '0.05'
        def mode = context.params.run_mode ?: 'default'

        // Pass as channel value
        process_data(Channel.value(threshold))
}
```

**CLI Usage:**
```bash
# Override parameter
bv run ./project --param.threshold 0.01
```

---

### 4. Bool

**Description:** Boolean true/false flag

**YAML Syntax:**
```yaml
parameters:
  - name: enable_filtering
    type: Bool
    default: true
    description: Enable quality filtering
```

**Nextflow Conversion:**
```groovy
workflow USER {
    take:
        context

    main:
        def enableFiltering = context.params.enable_filtering ?: true

        if (enableFiltering) {
            // Run filtering step
        }
}
```

**CLI Usage:**
```bash
# Set boolean parameter
bv run ./project --param.enable_filtering true
bv run ./project --param.enable_filtering false
```

---

### 5. Optional Types (Type?)

**Description:** Any type can be optional by adding `?`

**YAML Syntax:**
```yaml
inputs:
  - name: required_data
    type: File
    description: Required input file

  - name: optional_config
    type: File?
    description: Optional configuration file

  - name: optional_reference
    type: Directory?
```

**Pipeline Binding:**
```yaml
# Optional can be omitted
with:
  required_data: inputs.data
  # optional_config not provided - OK

# Or provided
with:
  required_data: inputs.data
  optional_config: File(config.yaml)
```

**Nextflow Conversion:**
```groovy
workflow USER {
    take:
        context
        required_data
        optional_config   // May be null/empty

    main:
        // Check if optional input provided
        if (optional_config) {
            process_with_config(required_data, optional_config)
        } else {
            process_without_config(required_data)
        }
}
```

**CLI Usage:**
```bash
# Without optional
bv run ./project --required_data data.csv

# With optional
bv run ./project \
  --required_data data.csv \
  --optional_config config.yaml
```

---

## Type Constructors

### File Constructor

**Syntax:** `File(path)`

**Usage:**
```yaml
# Pipeline literal binding
with:
  reference: File(references/hg38.fa)
  config: File(/absolute/path/config.json)

# Output declaration
publish:
  results: File(analysis_results.csv)
```

**Examples:**
```yaml
# Relative path
input: File(data/sample.csv)

# Absolute path
input: File(/var/data/reference.fa)

# Filename only (in publish)
output: File(results.txt)
```

### Directory Constructor

**Syntax:** `Directory(path)`

**Usage:**
```yaml
# Pipeline literal binding
with:
  data_dir: Directory(./data)
  reference_dir: Directory(/mnt/references)
```

---

## Format Specifications

For `File` types, specify format for validation and documentation:

**Supported Formats:**
- `csv` - Comma-separated values
- `tsv` - Tab-separated values
- `json` - JSON data
- `txt` - Plain text
- `vcf` - Variant Call Format
- `bam` - Binary Alignment Map
- `fasta` / `fa` - FASTA sequences
- `fastq` / `fq` - FASTQ sequences

**YAML Syntax:**
```yaml
outputs:
  - name: results
    type: File
    format: csv        # Format specification
    path: results.csv
```

**Benefits:**
- Documentation (clarifies expected format)
- Validation (future: schema checking)
- Format detection (for SQL stores, etc.)

---

## Type Compatibility

### Compatible Types

```yaml
# Same type
File → File              ✅
Directory → Directory    ✅
String → String          ✅

# Optional compatibility
File → File?             ✅
File? → File             ✅
Directory → Directory?   ✅
Directory? → Directory   ✅

# Case insensitive
File → file              ✅
Directory → DIRECTORY    ✅
```

### Incompatible Types

```yaml
# Different types
File → Directory         ❌
String → File            ❌
Bool → String            ❌

# Type errors caught during validation
bv pipeline validate pipeline.yaml
# Error: Type mismatch - expected File but found Directory
```

---

## Project YAML Syntax

### Full Project Specification

```yaml
name: project-name                # Required
author: email@example.com         # Required
workflow: workflow.nf             # Required
template: dynamic-nextflow        # Required
version: 1.0.0                    # Required

assets:                           # Optional
  - script.py
  - schema.yaml

parameters:                       # Optional
  - name: threshold
    type: String
    default: "30"
    description: Quality threshold

  - name: enable_qc
    type: Bool
    default: true

inputs:                           # Required
  - name: samplesheet
    type: File
    format: csv
    description: Input sample sheet

  - name: data_dir
    type: Directory
    description: Data directory

  - name: optional_config
    type: File?
    description: Optional configuration

outputs:                          # Required
  - name: results
    type: File
    format: csv
    path: analysis_results.csv
    description: Analysis output

  - name: summary
    type: File
    path: summary.txt
```

---

## Pipeline YAML Syntax

### Full Pipeline Specification

```yaml
name: pipeline-name               # Required

inputs:                           # Optional - pipeline-level inputs
  input_name: Type
  optional_input: Type?
  detailed_input:
    type: Type
    default: Type(path)

context:                          # Optional - shared context
  literal:
    key: value
    setting: option
  from_json: config.json

steps:                            # Required
  - id: step_id                   # Required - unique step name
    uses: project-reference       # Required - project path or name

    with:                         # Required - input bindings
      input1: inputs.name
      input2: step.other.outputs.data
      input3: File(literal.csv)

    publish:                      # Optional - output aliases
      alias_name: File(path)
      output_name: Type(path)

    store:                        # Optional - database storage
      store_name:
        kind: sql
        destination: SQL()
        source: output_name
        table_name: table_{run_id}
        key_column: id_column
        overwrite: true
        format: csv
```

---

## Nextflow Integration

### Type Mapping Table

| YAML Type | Nextflow Input | Nextflow Processing | Example |
|-----------|----------------|---------------------|---------|
| `File` | `path` | Staged file object | `path input.csv` |
| `Directory` | `val` | String path | `val "/data/dir"` |
| `String` | `val` | String value | `val "threshold"` |
| `Bool` | `val` | Boolean value | `val true` |
| `File?` | `path` (optional) | May be empty | `path config` |

### Workflow Signature

**YAML:**
```yaml
inputs:
  - name: samplesheet
    type: File
  - name: data_dir
    type: Directory
```

**Nextflow:**
```groovy
workflow USER {
    take:
        context       // Always first - BiovaultContext
        samplesheet   // File - staged path object
        data_dir      // Directory - string path value
    // ...
}
```

### Process Integration

**Using File:**
```groovy
process analyze_samples {
    input:
        path sample_csv    // File type → path

    output:
        path 'results.csv'

    script:
    """
    analyze.py ${sample_csv} > results.csv
    """
}
```

**Using Directory:**
```groovy
process scan_directory {
    input:
        val data_dir       // Directory type → val

    output:
        path 'file_list.txt'

    script:
    """
    ls ${data_dir}/*.csv > file_list.txt
    """
}
```

**Using Parameters:**
```groovy
workflow USER {
    take:
        context
        input_data

    main:
        // Access parameters from context
        def threshold = context.params.threshold
        def enableQC = context.params.enable_qc

        // Use in processes
        process_data(input_data, Channel.value(threshold))
}
```

---

## Complete Examples

### Example 1: Simple File Processing

**project.yaml:**
```yaml
name: file-processor
author: user@example.com
workflow: workflow.nf
template: dynamic-nextflow
version: 1.0.0

parameters:
  - name: min_size
    type: String
    default: "100"

inputs:
  - name: input_file
    type: File
    format: csv

outputs:
  - name: processed_file
    type: File
    format: csv
    path: output.csv
```

**workflow.nf:**
```groovy
nextflow.enable.dsl=2

workflow USER {
    take:
        context
        input_file     // File → path object

    main:
        def minSize = context.params.min_size
        result_ch = process_file(input_file, Channel.value(minSize))

    emit:
        processed_file = result_ch
}

process process_file {
    publishDir params.results_dir, mode: 'copy'

    input:
        path csv_file
        val min_size

    output:
        path 'output.csv'

    script:
    """
    filter.py --input ${csv_file} --min-size ${min_size} --output output.csv
    """
}
```

### Example 2: Multi-Step Pipeline

**pipeline.yaml:**
```yaml
name: analysis-pipeline

inputs:
  raw_data: File
  reference: File?

steps:
  - id: preprocess
    uses: ./preprocessor
    with:
      input: inputs.raw_data
      reference: inputs.reference
    publish:
      cleaned_data: File(cleaned.csv)

  - id: analyze
    uses: ./analyzer
    with:
      data: step.preprocess.outputs.cleaned_data
    publish:
      results: File(analysis.csv)
    store:
      analysis_results:
        kind: sql
        source: results
        table_name: analysis_{run_id}
        key_column: sample_id
```

**Running:**
```bash
# Without optional reference
bv run pipeline.yaml \
  --set inputs.raw_data=samples.csv

# With optional reference
bv run pipeline.yaml \
  --set inputs.raw_data=samples.csv \
  --set inputs.reference=ref.fa
```

---

## Best Practices

### 1. Use Specific Types

```yaml
# Good
- name: samplesheet
  type: File
  format: csv

# Avoid
- name: data
  type: File  # What kind of file?
```

### 2. Document with Descriptions

```yaml
# Good
- name: quality_threshold
  type: String
  default: "30"
  description: Minimum phred quality score (0-40)

# Avoid
- name: threshold
  type: String
  default: "30"
```

### 3. Use Formats for Files

```yaml
# Good
outputs:
  - name: results
    type: File
    format: vcf
    path: variants.vcf

# Avoid
outputs:
  - name: results
    type: File
    path: output.txt  # What format?
```

### 4. Mark Optional When Appropriate

```yaml
# Good
inputs:
  - name: required_samples
    type: File

  - name: optional_reference
    type: File?

# Avoid - making everything optional
inputs:
  - name: samples
    type: File?  # Should this really be optional?
```

### 5. Use Type Constructors in Pipelines

```yaml
# Good
with:
  config: File(configs/default.yaml)
  data_dir: Directory(./data)

# Avoid - unclear type
with:
  config: configs/default.yaml  # Is this a string or file?
```

---

## Common Patterns

### Pattern 1: Forward Directory

```yaml
# Project 1 outputs
outputs:
  - name: filtered_data
    type: File
  - name: data_dir
    type: Directory

# Project 2 uses same directory
inputs:
  - name: filtered_data
    type: File
  - name: data_dir
    type: Directory

# Pipeline chains them
steps:
  - id: filter
    outputs:
      filtered_data: File(filtered.csv)
      data_dir: Directory

  - id: analyze
    with:
      filtered_data: step.filter.outputs.filtered_data
      data_dir: step.filter.outputs.data_dir  # Forward directory
```

### Pattern 2: Optional Configuration

```yaml
inputs:
  - name: required_data
    type: File
  - name: optional_config
    type: File?

# Nextflow handles gracefully
workflow USER {
    take:
        context
        required_data
        optional_config

    main:
        if (optional_config) {
            config_ch = Channel.fromPath(optional_config)
        } else {
            config_ch = Channel.value([:])  // Empty config
        }
}
```

### Pattern 3: Parameter Defaults

```yaml
parameters:
  - name: threshold
    type: String
    default: "30"

# Override in pipeline or CLI
bv run project --param.threshold 50
```

---

## Troubleshooting

### Type Errors

**Error:** `Type mismatch: expected File but found Directory`
```yaml
# Problem
with:
  input_file: step.previous.outputs.data_dir  # Directory

# Fix
with:
  input_file: step.previous.outputs.result_file  # File
```

**Error:** `Unknown type 'List'`
```yaml
# Problem
- name: items
  type: List  # Not supported yet

# Fix - use File with CSV/JSON
- name: items
  type: File
  format: csv
```

**Error:** `Required input not provided`
```bash
# Problem
bv run project  # Missing required input

# Fix
bv run project --input_file data.csv
```

---

## Summary

### Quick Reference

| Type | YAML | Nextflow | Filesystem |
|------|------|----------|------------|
| File | `type: File` | `path` | Single file |
| Directory | `type: Directory` | `val` | Folder path |
| String | `type: String` | `val` | Text value |
| Bool | `type: Bool` | `val` | true/false |
| Optional | `type: Type?` | Optional | May be absent |

### Type Constructor Syntax

```yaml
File(path)              # File literal
Directory(path)         # Directory literal
Type(value)             # Generic constructor
```

### Binding Syntax

```yaml
inputs.name                              # Pipeline input
step.step_id.outputs.output_name        # Step output
File(literal/path)                      # Literal value
```

---

## Next Steps

- See [Project System Guide](project-system-guide.md) for project creation
- See [Pipeline System Guide](pipeline-system-guide.md) for pipeline composition
- Check [Nextflow Design Docs](nextflow_project_pipeline_design.md) for advanced patterns
