# Sheet Template

The Sheet template is designed for processing tabular data from CSV or TSV files with flexible schema validation and column mapping.

## Overview

The Sheet template enables you to:
- Process multiple samples from a CSV/TSV samplesheet
- Define schemas for validation and type conversion
- Rename columns to canonical names
- Set default values for missing fields
- Automatically cast data types

## Usage

### Basic Usage

```bash
# Sheet template REQUIRES a CSV/TSV file path
bv run ./my-project /path/to/samplesheet.csv

# Using a TSV file
bv run ./my-project data.tsv

# With options
bv run ./my-project participants.csv --dry-run --with-docker
```

**Note:** Unlike other templates, sheet templates require you to explicitly provide the samplesheet path. There is no default - this ensures you're always explicit about which data you're processing.

## Project Structure

A sheet-based project requires:

```
my-project/
├── project.yaml          # Must specify template: sheet
├── workflow.nf           # Your workflow implementation
└── assets/
    └── schema.yaml       # Column validation rules (optional but recommended)
```

**Note:** The samplesheet (CSV/TSV file) is provided at runtime, not stored in the project.

### project.yaml

```yaml
name: my-sheet-project
author: user@example.com
workflow: workflow.nf
template: sheet
assets:
  - schema.yaml  # Include schema.yaml in assets for hashing/tracking
```

### assets/schema.yaml

The schema file defines validation rules and transformations:

```yaml
# Required columns (will error if missing)
required:
  - participant_id
  - genotype_file_path

# Rename columns (from CSV header to canonical name)
rename:
  geno_path: genotype_file_path
  sample: participant_id

# Default values for missing/empty fields
defaults:
  age: 40
  ref_version: GRCh38

# Type casting (string, int, float, bool, path)
types:
  participant_id: string
  genotype_file_path: path
  weight: float
  height: float
  age: int
  is_case: bool
```

### workflow.nf

Your workflow receives:
- `rows_ch`: Channel of row maps (each row as a Map)
- `mapping_ch`: Single map with headers/types/schema info
- `assets_dir_ch`: Path to assets directory
- `results_dir`: Results directory path

Example workflow:

```nextflow
nextflow.enable.dsl=2

workflow USER {
    take:
    rows_ch         // Channel of row maps
    mapping_ch      // Schema metadata
    assets_dir_ch   // Assets directory
    results_dir     // Results directory

    main:
    // Process each row
    rows_ch
        .map { row ->
            // Access fields by name
            def id = row.participant_id
            def geno = row.genotype_file_path
            // Return tuple for processing
            tuple(id, geno)
        }
        .set { samples_ch }

    // Run your analysis
    process_samples(samples_ch)

    emit:
    done = process_samples.out
}

process process_samples {
    input:
    tuple val(id), path(genotype)

    output:
    path "${id}_results.txt"

    script:
    """
    echo "Processing ${id}" > ${id}_results.txt
    # Your analysis here
    """
}
```

## Samplesheet Format

### CSV Format
```csv
participant_id,genotype_file_path,age,weight
SAMPLE001,/data/sample001.vcf,35,70.5
SAMPLE002,/data/sample002.vcf,42,65.0
SAMPLE003,/data/sample003.vcf,,80.2
```

### TSV Format
Use `.tsv` extension for tab-separated values:
```tsv
participant_id	genotype_file_path	age	weight
SAMPLE001	/data/sample001.vcf	35	70.5
SAMPLE002	/data/sample002.vcf	42	65.0
```

## Type System

The sheet template supports automatic type conversion:

- **string**: Default type, no conversion
- **int**: Converts to integer
- **float**: Converts to decimal number
- **bool**: Converts "true", "1", "yes", "y" to true; others to false
- **path**: Converts to Nextflow file/path object

## Schema Processing Order

1. **Rename**: Column headers are renamed according to `rename` mapping
2. **Defaults**: Empty/missing values are filled with defaults
3. **Required**: Check that all required fields have values
4. **Types**: Cast values to specified types

## Examples

### Example 1: Simple Genotype Processing

**participants.csv:**
```csv
sample,geno_file
S001,/data/s001.vcf
S002,/data/s002.vcf
```

**schema.yaml:**
```yaml
rename:
  sample: participant_id
  geno_file: genotype_path
required:
  - participant_id
  - genotype_path
types:
  genotype_path: path
```

### Example 2: Clinical Data with Defaults

**participants.csv:**
```csv
id,age,is_case,bmi
P001,45,true,25.3
P002,,false,28.1
P003,52,,22.0
```

**schema.yaml:**
```yaml
rename:
  id: participant_id
defaults:
  age: 50
  is_case: false
types:
  age: int
  is_case: bool
  bmi: float
```

## Error Handling

The template will fail if:
- Required columns are missing from the samplesheet
- Required fields have empty values (after defaults)
- Type conversion fails (e.g., non-numeric value for int type)
- The samplesheet file is empty or unreadable

## Tips

1. **File Paths**: Use `path` type for file inputs to enable proper staging
2. **Validation**: Use `required` to ensure critical fields are present
3. **Defaults**: Set sensible defaults to handle missing optional data
4. **Column Names**: Use `rename` to map user-friendly names to canonical ones
5. **Extension**: Use `.tsv` for tab-separated files, `.csv` for comma-separated

## See Also

- [sheet-add example](../examples/sheet-add/) - Example project using sheet template
- [CSV Batch Processing](csv_batch_processing.md) - For processing multiple participants with standard templates