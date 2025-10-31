# BioVault SQL System Guide

## Overview

BioVault provides an integrated SQL system for storing, querying, and managing tabular data. The system consists of:

1. **SQLite Database** - Built-in database for operational and result data
2. **Pipeline SQL Stores** - Automatic storage of pipeline outputs to database tables
3. **SQL CLI Commands** - Tools for querying, importing, and managing data

All data is stored in a local SQLite database at `{BioVault home}/biovault.db` (default: `~/Desktop/BioVault/biovault.db`, or `BIOVAULT_HOME/biovault.db` if set).

---

## Database Structure

### Core Tables

BioVault maintains core operational tables:

| Table | Purpose | Protected |
|-------|---------|-----------|
| `participants` | Participant registry | ✅ Yes |
| `projects` | Registered projects | ✅ Yes |
| `messages` | Inter-datasite messages | ✅ Yes |
| `message_threads` | Message threading | ✅ Yes |

**Protected tables cannot be overwritten or dropped.**

### Result Tables

Pipeline results are stored in tables with the `z_results_` prefix:

| Table Pattern | Source | Protected |
|---------------|--------|-----------|
| `z_results_*` | Pipeline SQL stores | ❌ No (can be dropped) |
| `user_*` | Manual imports (future) | ❌ No (can be dropped) |

**Benefits of `z_results_` prefix:**
- Tables sort last alphabetically
- Easy to identify pipeline outputs
- Simple cleanup of old runs

---

## Pipeline SQL Stores

### Overview

Pipeline SQL stores automatically save step outputs to the database after execution.

**Use cases:**
- Persist analysis results for querying
- Track metrics across pipeline runs
- Enable SQL-based downstream analysis
- Archive results without keeping large files

### Configuration

Add a `store:` section to any pipeline step:

```yaml
steps:
  - id: analyze
    uses: ./analysis-project
    with:
      data: inputs.samplesheet
    publish:
      results: File(analysis.csv)
    store:
      analysis_results:           # Store identifier
        kind: sql                 # Store type (currently only 'sql')
        destination: SQL()        # Target database
        source: results           # Output to store
        table_name: analysis_{run_id}
        key_column: participant_id
        overwrite: true           # Optional, default: true
        format: csv               # Optional, auto-detected
```

### Store Fields

#### Required Fields

**`kind`** - Store type
- Currently only `sql` is supported
- Future: `s3`, `gcs`, `api`, etc.

**`source`** - Output name to store
- Must match a name in `publish:` or project outputs
- References the CSV/TSV file to import

#### Optional Fields

**`destination`** - Database target
- `SQL()` or omit → BioVault SQLite database
- `SQL(url:...)` → External database (not yet supported)

**`table_name`** - Target table name
- Supports `{run_id}` substitution for unique runs
- Defaults to `{store_name}_{run_id}`
- Automatically prefixed with `z_results_`

Example:
```yaml
table_name: analysis_{run_id}
# Creates: z_results_analysis_20251023120000
```

**`key_column`** - Primary key column
- Optional - specifies which column is the primary key
- Column must exist in the source CSV/TSV
- Case-insensitive matching

**`overwrite`** - Drop existing table
- Default: `true` (drop and recreate)
- `false` - fail if table exists

**`format`** - Data format
- `csv` or `tsv`
- Auto-detected from file extension if omitted

### Example Pipeline

**pipeline.yaml:**
```yaml
name: variant-analysis

inputs:
  samplesheet: File
  data_dir: Directory

steps:
  - id: count_variants
    uses: ./variant-counter
    with:
      samplesheet: inputs.samplesheet
      data_dir: inputs.data_dir
    publish:
      variant_counts: File(counts.csv)
    store:
      variant_db:
        kind: sql
        destination: SQL()
        source: variant_counts
        table_name: variants_{run_id}
        key_column: sample_id

  - id: quality_stats
    uses: ./quality-analyzer
    with:
      samplesheet: step.count_variants.outputs.variant_counts
    publish:
      quality_metrics: File(metrics.csv)
    store:
      quality_db:
        kind: sql
        source: quality_metrics
        table_name: quality_metrics_{run_id}
```

**Running:**
```bash
bv run pipeline.yaml \
  --set inputs.samplesheet=samples.csv \
  --set inputs.data_dir=./data
```

**Output:**
```
💾 Stored 'variant_db' output 'variant_counts' into table z_results_variants_20251023120000 (rows: 150).
    source: /path/to/results/count_variants/counts.csv
    database: /Users/user/Desktop/BioVault/biovault.db

💾 Stored 'quality_db' output 'quality_metrics' into table z_results_quality_metrics_20251023120000 (rows: 150).
    source: /path/to/results/quality_stats/metrics.csv
    database: /Users/user/Desktop/BioVault/biovault.db
```

---

## SQL CLI Commands

### Available Commands

```bash
bv sql <command> [options]
```

| Command | Purpose | Status |
|---------|---------|--------|
| `tables` | List all tables | ✅ Implemented |
| `structure` | Show table schema | ✅ Implemented |
| `run` | Execute SQL query | ✅ Implemented (read-only) |
| `import` | Import CSV/TSV to table | 🚧 Planned |
| `drop` | Delete user tables | 🚧 Planned |

---

### 1. List Tables

View all tables in the database:

```bash
bv sql tables
```

**Example Output:**
```
participants
projects
messages
message_threads
z_results_analysis_20251023120000
z_results_variants_20251023115500
```

**Use cases:**
- Quick reference of available data
- Find pipeline result tables
- List tables for AI/scripts

---

### 2. Table Structure

View schema information for tables:

```bash
# Single table
bv sql structure <table_name>

# All tables
bv sql structure
bv sql structure all
```

**Example Output:**
```
TABLE: z_results_variants_20251023120000
├─ sample_id: TEXT PRIMARY KEY
├─ variant_count: TEXT
├─ quality_score: TEXT
└─ created_at: TEXT
```

**JSON Output:**
```bash
bv sql structure z_results_variants_20251023120000 --json
```

```json
{
  "table": "z_results_variants_20251023120000",
  "columns": [
    {"name": "sample_id", "type": "TEXT", "constraints": ["PRIMARY KEY"]},
    {"name": "variant_count", "type": "TEXT"},
    {"name": "quality_score", "type": "TEXT"}
  ]
}
```

---

### 3. Run Queries

Execute SQL queries on the database:

```bash
bv sql run "<query>" [OPTIONS]
```

**Options:**
- `--output <path>` or `-o <path>` - Export results to file
- `--format <csv|tsv|json>` - Output format (default: csv)
- `--plain` - Plain output without formatting

**Basic Query:**
```bash
# Display results
bv sql run "SELECT * FROM z_results_variants_20251023120000 LIMIT 10"
```

**Example Output:**
```
sample_id    | variant_count | quality_score
-------------|---------------|---------------
SAMPLE_001   | 42            | 98.5
SAMPLE_002   | 38            | 97.2
SAMPLE_003   | 51            | 99.1
...
(10 rows)
```

**Export to CSV:**
```bash
bv sql run "SELECT * FROM z_results_variants_20251023120000" \
  --output results.csv
```

**Export as TSV:**
```bash
bv sql run "SELECT sample_id, variant_count FROM z_results_variants_20251023120000" \
  --output data.tsv \
  --format tsv
```

**JSON Output:**
```bash
bv sql run "SELECT * FROM z_results_variants_20251023120000" \
  --format json
```

---

## Data Import & Storage Process

### How Pipeline SQL Stores Work

When a pipeline step completes:

1. **Output Detection** - CLI identifies outputs marked for storage
2. **File Reading** - Reads CSV/TSV from step results directory
3. **Table Creation** - Creates table with sanitized column names
4. **Schema Generation** - Auto-generates schema from CSV headers
5. **Primary Key** - Sets PRIMARY KEY on `key_column` if specified
6. **Data Import** - Inserts all rows using transaction
7. **Table Prefixing** - Final table name prefixed with `z_results_`

### CSV Processing

**Header Handling:**
```csv
participant_id,genotype_file_path,line_count
PID00001,case_0001_canonical.txt,4
PID00002,case_0002_canonical.txt,4
```

**Sanitization Rules:**
- Uppercase → lowercase: `ParticipantID` → `participantid`
- Spaces → underscores: `Sample ID` → `sample_id`
- Special chars → underscores: `Sample-ID@Test` → `sample_id_test`
- Multiple underscores → single: `sample___id` → `sample_id`
- Leading/trailing underscores trimmed
- Starts with number → prefix `t_`: `123abc` → `t_123abc`
- Empty column names → `col1`, `col2`, etc.

**Example:**
```csv
Participant ID,Genotype-File,Line Count!
```
Becomes:
```sql
CREATE TABLE z_results_table_name (
  "participant_id" TEXT PRIMARY KEY,
  "genotype_file" TEXT,
  "line_count" TEXT
)
```

### Type Handling

**Current:** All columns stored as `TEXT`
- Preserves data exactly as provided
- No type conversion errors
- Flexible for mixed content

**Future:** Type inference
- Detect numeric columns → `INTEGER`, `REAL`
- Detect dates → `TIMESTAMP`
- Boolean detection → `INTEGER` (0/1)

---

## Querying Pipeline Results

### Common Patterns

**List All Result Tables:**
```bash
bv sql run "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'z_results_%'"
```

**Recent Runs:**
```bash
# Table names include timestamp: z_results_analysis_20251023120000
bv sql run "SELECT name FROM sqlite_master
  WHERE type='table' AND name LIKE 'z_results_%'
  ORDER BY name DESC LIMIT 10"
```

**Aggregate Across Runs:**
```bash
# Join multiple run tables (requires identical schemas)
bv sql run "
  SELECT '20251023120000' as run_id, * FROM z_results_analysis_20251023120000
  UNION ALL
  SELECT '20251022110000' as run_id, * FROM z_results_analysis_20251022110000
" --output combined_results.csv
```

**Filter Results:**
```bash
bv sql run "
  SELECT * FROM z_results_variants_20251023120000
  WHERE variant_count > 40
  AND quality_score > 95.0
"
```

**Summary Statistics:**
```bash
bv sql run "
  SELECT
    COUNT(*) as total_samples,
    AVG(CAST(variant_count AS REAL)) as avg_variants,
    MAX(CAST(quality_score AS REAL)) as max_quality
  FROM z_results_variants_20251023120000
"
```

---

## Table Management

### Cleanup Old Results

Pipeline results accumulate over time. Clean up manually:

```bash
# List old result tables
bv sql run "SELECT name FROM sqlite_master
  WHERE type='table' AND name LIKE 'z_results_%'
  ORDER BY name"

# Drop specific table (future)
bv sql drop z_results_old_analysis_20251020000000
```

**Current Workaround:**
```bash
# Direct SQL for now
bv sql run "DROP TABLE z_results_old_table_name"
```

### Archive Results

Export results before cleanup:

```bash
# Export to CSV
bv sql run "SELECT * FROM z_results_analysis_20251020000000" \
  --output archived_results_20251020.csv

# Then drop
bv sql run "DROP TABLE z_results_analysis_20251020000000"
```

---

## Best Practices

### 1. Meaningful Table Names

```yaml
# Good - descriptive with run_id
store:
  variant_analysis:
    table_name: variant_counts_{run_id}

# Avoid - static names
store:
  results:
    table_name: results  # Overwrites every run!
```

### 2. Use Primary Keys

```yaml
# Good - enables efficient lookups
store:
  analysis_db:
    key_column: participant_id

# OK - no primary key
store:
  logs_db:
    # No key_column - table has no primary key
```

### 3. Consistent Column Names

Ensure CSV outputs use consistent, SQL-friendly names:

```python
# Good
df.columns = ['participant_id', 'variant_count', 'quality_score']

# Avoid
df.columns = ['Participant ID', 'Variant-Count', 'Quality Score!']
```

### 4. Document Table Schemas

Add comments to clarify table contents:

```yaml
store:
  variant_counts:
    # Stores per-sample variant counts with quality metrics
    # Schema: sample_id (PK), variant_count (int), quality_score (float)
    kind: sql
    source: variant_results
    table_name: variant_counts_{run_id}
    key_column: sample_id
```

### 5. Regular Cleanup

Automate result cleanup:

```bash
#!/bin/bash
# cleanup_old_results.sh

# Keep last 30 days of results
CUTOFF_DATE=$(date -d '30 days ago' +%Y%m%d)

bv sql run "SELECT name FROM sqlite_master
  WHERE type='table' AND name LIKE 'z_results_%'" --plain | \
while read table; do
  # Extract date from table name (assumes format: z_results_name_YYYYMMDDHHMMSS)
  table_date=$(echo "$table" | grep -oE '[0-9]{14}' | cut -c1-8)

  if [[ "$table_date" < "$CUTOFF_DATE" ]]; then
    echo "Dropping old table: $table"
    bv sql run "DROP TABLE $table"
  fi
done
```

---

## Security & Safety

### SQL Injection Prevention

The system prevents SQL injection:

```bash
# ✅ Safe - parameterized queries (future)
bv sql run "SELECT * FROM users WHERE id = ?" --params "42"

# ⚠️ Current - validates query structure
bv sql run "SELECT * FROM users WHERE id = 42"

# ❌ Blocked - dangerous patterns detected
bv sql run "SELECT * FROM users WHERE id = 1; DROP TABLE users"
```

### Protected Tables

Core BioVault tables are protected:

```bash
# ❌ Error - cannot drop protected tables
bv sql run "DROP TABLE participants"
# Error: Table 'participants' is protected and cannot be dropped

# ✅ OK - can drop result tables
bv sql run "DROP TABLE z_results_old_analysis_20251020000000"
```

### Read-Only by Default

Currently all `bv sql run` commands are read-only (SELECT only).

**Future write operations:**
```bash
# Write operations will require flag
bv sql run "INSERT INTO user_data ..." --allow-write

# DDL operations require confirmation
bv sql run "DROP TABLE user_table" --allow-ddl
```

---

## Troubleshooting

### Table Not Found

```
Error: no such table: z_results_analysis
```

**Solution:** Check exact table name with timestamp:
```bash
bv sql tables | grep analysis
# Shows: z_results_analysis_20251023120000
```

### Column Name Issues

```
Error: no such column: Participant-ID
```

**Solution:** Columns are sanitized - use lowercase with underscores:
```bash
# Instead of:
bv sql run "SELECT 'Participant-ID' FROM ..."

# Use:
bv sql run "SELECT participant_id FROM ..."
```

### Primary Key Constraint Error

```
Error: UNIQUE constraint failed: z_results_table.sample_id
```

**Solution:** Duplicate values in key_column:
- Check source CSV for duplicate IDs
- Remove duplicates before storing
- Or omit `key_column` if duplicates are valid

### Empty Results

```
(0 rows)
```

**Solution:** Verify data was stored:
```bash
# Check table exists and has data
bv sql run "SELECT COUNT(*) FROM z_results_table_name"
```

---

## Future Enhancements

### Planned Features

**Advanced Querying:**
- Joins across pipeline runs
- Saved query templates
- Query builder UI
- Performance analysis

**Data Import:**
```bash
# Import custom CSV to user table
bv sql import data.csv --table user_custom_data

# With column mapping
bv sql import survey.csv \
  --table survey_results \
  --map "subject_id:participant_id,response:answer"

# Participant validation
bv sql import data.csv \
  --table results \
  --participant-col sample_id \
  --on-mismatch skip  # or 'create' or 'error'
```

**Backup & Migration:**
```bash
# Full database backup
bv sql backup --output biovault_backup.db

# Export all result tables
bv sql export-all --output-dir ./sql_exports
```

**External Databases:**
```yaml
# Store to PostgreSQL (future)
store:
  external_results:
    kind: sql
    destination: SQL(url:postgresql://localhost/mydb)
    source: results
```

---

## Complete Example

### Multi-Step Pipeline with SQL Storage

**pipeline.yaml:**
```yaml
name: comprehensive-analysis

inputs:
  raw_samples: File
  reference_data: Directory

steps:
  - id: quality_filter
    uses: ./qc-filter
    with:
      samples: inputs.raw_samples
      reference: inputs.reference_data
    publish:
      filtered_samples: File(qc_passed.csv)
    store:
      qc_results:
        kind: sql
        source: filtered_samples
        table_name: qc_results_{run_id}
        key_column: sample_id

  - id: variant_calling
    uses: ./variant-caller
    with:
      samples: step.quality_filter.outputs.filtered_samples
      reference: inputs.reference_data
    publish:
      variants: File(variants.csv)
    store:
      variant_results:
        kind: sql
        source: variants
        table_name: variants_{run_id}
        key_column: sample_id

  - id: annotation
    uses: ./annotator
    with:
      variants: step.variant_calling.outputs.variants
    publish:
      annotated_variants: File(annotated.csv)
    store:
      final_results:
        kind: sql
        source: annotated_variants
        table_name: annotated_variants_{run_id}
        key_column: variant_id
```

**Run Pipeline:**
```bash
bv run pipeline.yaml \
  --set inputs.raw_samples=samples.csv \
  --set inputs.reference_data=./ref_data \
  --results-dir ./analysis_results
```

**Query Results:**
```bash
# List created tables
bv sql tables | grep z_results

# View QC statistics
bv sql run "SELECT
    COUNT(*) as total_passed,
    AVG(CAST(quality_score AS REAL)) as avg_quality
  FROM z_results_qc_results_20251023120000"

# Export variants for downstream analysis
bv sql run "SELECT * FROM z_results_variants_20251023120000" \
  --output variants_for_analysis.csv

# Join QC and variant results
bv sql run "
  SELECT
    qc.sample_id,
    qc.quality_score,
    COUNT(v.variant_id) as variant_count
  FROM z_results_qc_results_20251023120000 qc
  LEFT JOIN z_results_variants_20251023120000 v
    ON qc.sample_id = v.sample_id
  GROUP BY qc.sample_id
" --output combined_metrics.csv
```

---

## Summary

### Key Concepts

✅ **Pipeline SQL Stores** - Automatic database storage of outputs
✅ **Table Prefixing** - `z_results_` prefix for easy identification
✅ **Run Isolation** - `{run_id}` in table names separates runs
✅ **SQL CLI** - Query and export tools
✅ **Protected Tables** - Core tables cannot be modified
✅ **Type Safety** - Automatic schema generation from CSV

### Quick Reference

```bash
# View tables
bv sql tables

# Show schema
bv sql structure <table_name>

# Query data
bv sql run "SELECT * FROM <table_name>"

# Export results
bv sql run "SELECT * FROM <table_name>" -o output.csv

# Pipeline with SQL store
steps:
  - id: analyze
    store:
      results_db:
        kind: sql
        source: output_name
        table_name: results_{run_id}
        key_column: id_column
```

---

## Next Steps

- See [Pipeline System Guide](pipeline-system-guide.md) for pipeline creation
- See [YAML Types Reference](yaml-types-reference.md) for type specifications
- Check [Project System Guide](project-system-guide.md) for project development
