# BioVault

BioVault is a free, open-source, permissionless network for collaborative genomics.

Built with end-to-end encryption, secure enclaves, and data visitation, BioVault lets researchers and participants share insights without ever sharing raw data.

https://biovault.net/

## Quick Install (One-liner)

```bash
curl -sSL https://raw.githubusercontent.com/openmined/biovault/main/install.sh | bash
```

## Prerequisites
- [SyftBox](https://syftbox.net)
- [NextFlow](https://www.nextflow.io)
  - [Java 17+](https://openjdk.java.net/)
- [Docker](https://www.docker.com) (optional)

## Setup
Run `bv check` and make sure you have the depenencies listed below.
```
bv check

BioVault Dependency Check
=========================

Checking java...  (version 23)✓ Found
Checking docker... ✓ Found (running)
Checking nextflow... ✓ Found
Checking syftbox... ✓ Found

=========================
✓ All dependencies satisfied!
```

## Automatic Setup
You can `bv setup` on some systems such as macOS and Google Colab and `bv` will help you to install the dependencies.

## SyftBox
SyftBox requires setup and authentication.

## Tutorials:
- [1) Hello World](tutorials/1_hello_world.md)
- [2) Submit Your Project](tutorials/2_submit_your_project.md)
- [3) Create a Biobank](tutorials/3_create_biobank.md)

## Documentation

- [Development Guide](DEV.md) - Setup and testing instructions
- [Security](SECURITY.md) - How BioVault protects your data with SyftBox permissions

## Development

For development setup and commands, see [DEV.md](DEV.md).



## CLI Overview

The `bv` CLI provides commands to manage BioVault projects, data, messaging, and utilities.

Global flags
- `-v, --verbose` Increase log verbosity
- `--config <path>` Use a specific config file

Top-level commands
- `bv update` Check for updates and install the latest
- `bv init [email]` Initialize a new BioVault repo; email is optional (detected from `SYFTBOX_EMAIL` if omitted)
- `bv info` Show system information
- `bv check` Check for required dependencies
- `bv setup` Setup environment for known systems (e.g., Google Colab)
- `bv project create [--name <name>] [--folder <path>]` Create a new project scaffold
- `bv run <project_folder> <participant_source> [--test] [--download] [--dry-run] [--with-docker=<bool>] [--work-dir <dir>] [--resume]`
  - `participant_source` can be a local file path, Syft URL, or HTTP URL (with optional `#fragment`)
  - `--with-docker` defaults to `true`
- `bv sample-data fetch [--participant-ids id1,id2,...] [--all]` Fetch sample data
- `bv sample-data list` List available sample data
- `bv participant add [--id <ID>] [--aligned <file>]` Add a participant record
- `bv participant list` List participants
- `bv participant delete <ID>` Delete a participant
- `bv participant validate [--id <ID>]` Validate participant files (all if omitted)
- `bv biobank list` List biobanks in SyftBox
- `bv biobank publish [--participant-id <ID>] [--all] [--http-relay-servers host1,host2,...]` Publish participants
- `bv biobank unpublish [--participant-id <ID>] [--all]` Unpublish participants
- `bv config email <email>` Set email address
- `bv config syftbox [--path <config.json>]` Set SyftBox config path
- `bv fastq combine <input_folder> <output_file> [--validate] [--no-prompt] [--stats-format tsv|yaml|json]` Combine/validate FASTQ files
- `bv submit <project_path> <destination>` Submit a project (destination is datasite email or full Syft URL)
- `bv samplesheet create <input_dir> <output_file> [--file_filter <pattern>] [--extract_cols <pattern>] [--ignore]` Create sample sheet CSV from files

Inbox and messaging
- `bv inbox` Interactive inbox (default; uses single-key shortcuts)
  - Shortcuts: `?`/`h` Help, `n` New, `s` Sync, `v` Change view, `q` Quit, `1..5` Tabs (Inbox, Sent, All, Unread, Projects)
  - Arrow keys navigate; Enter opens the selected message or Quit
- `bv inbox --plain [--sent] [--all] [--unread] [--projects] [--type <text|project|request>] [--from <sender>] [--search <term>]`
  - Non-interactive list output with filters
- `bv message send <recipient> <message> [-s|--subject <subject>]` Send a message
- `bv message reply <message_id> <body>` Reply to a message
- `bv message read <message_id>` Read a specific message
- `bv message delete <message_id>` Delete a message
- `bv message list [--unread]` List messages (optionally only unread)
- `bv message thread <thread_id>` View a message thread
- `bv message sync` Sync messages (check for new and update ACKs)

Examples
- Initialize and set email: `bv init you@example.com`
- Create a new project: `bv project create --name demo --folder ./demo`
- Run a project with test data: `bv run ./demo participants.yaml --test --download`
- Combine FASTQs: `bv fastq combine ./fastq_pass ./combined/output.fastq.gz --validate`
- Interactive inbox: `bv inbox` (press `?` for shortcuts)
- Plain inbox list: `bv inbox --plain --unread`
- Create sample sheet from genotype files:
  ```bash
  # Extract participant IDs from filenames matching a pattern
  bv samplesheet create test_dir output.csv --extract_cols="{participant_id}_X_X_GSAv3-DTC_GRCh38-{date}.txt"

  # Example with files: 103704_X_X_GSAv3-DTC_GRCh38-07-01-2025.txt
  # Produces CSV:
  # participant_id,genotype_file_path
  # 103704,/absolute/path/test_dir/103704_X_X_GSAv3-DTC_GRCh38-07-01-2025.txt
  ```

## File Import Workflow

The `bv files` commands provide a flexible workflow for importing genomic data files with automatic participant ID extraction and file type detection.

### Complete Import Example

#### 1. Scan directory to see what file types are available

```bash
bv files scan /path/to/data
```

Output:
```
📊 Scan Results: /path/to/data

Extensions Found:
  .txt  323 files    6701.8 MB
  .csv  4 files    32.6 MB

Total: 332 files
```

#### 2. Suggest patterns for extracting participant IDs from file paths

```bash
bv files suggest-patterns /path/to/data --ext .txt
```

Output:
```
🔍 Detected Patterns:

1. {parent} - Directory name as participant ID
   Example: huE922FC/...
   Sample extractions:
     huE922FC/... → participant ID: huE922FC
     huBF0F93/... → participant ID: huBF0F93
```

#### 3. Preview import with dry-run

```bash
bv files import /path/to/data --ext .txt --pattern {parent} --dry-run
```

Output shows sample participant ID extractions without importing.

#### 4. Export file list to CSV with pattern-based participant ID extraction

```bash
bv files export-csv /path/to/data --ext .txt --pattern {parent} -o genotype-files.csv
```

Output:
```
📊 Found 323 files
✓ Exported 323 files to genotype-files.csv
```

#### 5. Detect file types and update CSV

```bash
bv files detect-csv genotype-files.csv -o genotype-files.csv
```

Output:
```
🔍 Detecting file types from genotype-files.csv
📋 Processing 323 files
🔍 Detecting... 323/323
✓ Updated CSV written to genotype-files.csv
```

#### 6. Import files using the CSV

```bash
bv files import-csv genotype-files.csv
```

Output:
```
📋 CSV Import Preview: genotype-files.csv
  Files to import: 323
```

### Pattern Examples

- `{parent}` - Use parent directory name as participant ID
- `{filename}` - Use filename as participant ID
- Custom patterns can extract from any part of the file path

### Available Commands

- `bv files scan <path>` - Scan directory and show file type statistics
- `bv files suggest-patterns <path> --ext <extension>` - Analyze files and suggest participant ID extraction patterns
- `bv files import <path> --ext <ext> --pattern <pattern> [--dry-run]` - Preview or import files with pattern
- `bv files export-csv <path> --ext <ext> --pattern <pattern> -o <output.csv>` - Export file list with participant IDs to CSV
- `bv files detect-csv <input.csv> -o <output.csv>` - Detect file types and update CSV
- `bv files import-csv <file.csv>` - Import files from CSV into BioVault database

## SyftBox VirtualEnv
If you need to run multiple syftbox instances checkout `sbenv` which will help you to isolate them on your machine:
https://github.com/openmined/sbenv

BioVault can auto detect when its in an `sbenv activate` environment and will target that isolated syftbox for all its usage.
