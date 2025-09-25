# CSV Batch Processing

The BioVault CLI now supports batch processing of multiple participants from a CSV file. This feature allows you to run the same workflow on multiple samples with different parameters in a single command.

## Basic Usage

```bash
bv run <project_folder> <csv_file>
```

Example:
```bash
bv run ./my-project participants.csv
```

## CSV File Formats

### 1. Auto-mapped CSV (Column names match parameter names)

When your CSV column names exactly match the workflow parameter names:

```csv
participant_id,ref_version,ref,ref_index,aligned,aligned_index
NA12878,GRCh38,/data/reference/GRCh38.fa,/data/reference/GRCh38.fa.fai,/data/samples/NA12878.cram,/data/samples/NA12878.cram.crai
NA12891,GRCh38,/data/reference/GRCh38.fa,/data/reference/GRCh38.fa.fai,/data/samples/NA12891.cram,/data/samples/NA12891.cram.crai
```

### 2. Custom Column Names with Mapping File

When your CSV has custom column names, provide a mapping file:

```bash
bv run ./my-project custom_data.csv:mapping.yaml
```

**CSV file (custom_data.csv):**
```csv
sample_name,genome_build,reference_genome,alignment_file
SAMPLE001,hg38,/lab/ref/hg38.fa,/lab/alignments/SAMPLE001.cram
SAMPLE002,hg38,/lab/ref/hg38.fa,/lab/alignments/SAMPLE002.cram
```

**Mapping file (mapping.yaml):**
```yaml
mappings:
  sample_name: participant_id
  genome_build: ref_version
  reference_genome: ref
  alignment_file: aligned
```

### 3. Inline Mapping (Column headers with colon syntax)

Use `column_name:parameter_name` format in the CSV header:

```csv
sample:participant_id,build:ref_version,ref_file:ref,cram_file:aligned
TestSample1,GRCh38,/data/GRCh38.fa,/data/test1.cram
TestSample2,GRCh38,/data/GRCh38.fa,/data/test2.cram
```

## Supported Parameters

The following standard parameters are automatically mapped:

- `participant_id` or `id` - Sample/participant identifier
- `ref_version` - Reference genome version (e.g., GRCh38, GRCh37)
- `ref` or `ref_path` - Path to reference genome FASTA file
- `ref_index` - Path to reference genome index file
- `aligned` - Path to aligned reads file (BAM/CRAM)
- `aligned_index` - Path to alignment index file
- `snp` - Path to SNP/variant file (VCF, 23andMe, etc.)

## Custom Parameters

Any additional columns in your CSV can be passed as custom parameters to the workflow. These will be available in your workflow.nf file as `params.<parameter_name>`.

Example CSV with custom parameters:
```csv
participant_id,ref_version,aligned,phenotype,batch_id
SAMPLE001,GRCh38,/data/s001.cram,case,batch1
SAMPLE002,GRCh38,/data/s002.cram,control,batch1
```

The `phenotype` and `batch_id` columns will be passed to the workflow as `params.phenotype` and `params.batch_id`.

## Command Options

All standard `bv run` options work with CSV batch processing:

```bash
# Dry run to preview commands
bv run ./my-project participants.csv --dry-run

# Use Docker for execution
bv run ./my-project participants.csv --with-docker

# Specify custom template
bv run ./my-project participants.csv --template custom_template

# Resume failed runs
bv run ./my-project participants.csv --resume

# Auto-download files
bv run ./my-project participants.csv --download
```

## Examples

See the `examples/` directory for complete examples:

- `batch_participants.csv` - Standard CRAM processing batch
- `custom_batch.csv` + `custom_mapping.yaml` - Custom column names with mapping
- `inline_mapping.csv` - Inline column mapping
- `snp_batch.csv` - SNP data batch processing

## Processing Order

- Each row in the CSV is processed sequentially
- The workflow runs completely for one participant before moving to the next
- Results are stored in separate directories for each participant
- If a run fails, subsequent rows will still be processed

## Tips

1. **Test First**: Use `--dry-run` to preview the commands that will be executed
2. **Consistent Paths**: Ensure all file paths in the CSV are accessible from where you run the command
3. **Template Compatibility**: Make sure your workflow template supports the parameters in your CSV
4. **Error Handling**: Check the results directory for each participant to identify any failures
5. **Large Batches**: For very large batches, consider splitting the CSV and running multiple instances

## Troubleshooting

**Issue: "CSV file not found"**
- Ensure the CSV file path is correct and the file exists

**Issue: "Unknown parameter"**
- Check that custom parameters are defined in your workflow.nf file

**Issue: "File not found" during execution**
- Verify all file paths in the CSV are absolute paths or relative to the working directory

**Issue: Mapping not working**
- Ensure mapping file is valid YAML with a `mappings:` section
- Check column names match exactly (case-sensitive)