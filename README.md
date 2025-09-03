# BioVault

## Quick Install (One-liner)

```bash
curl -sSL https://raw.githubusercontent.com/openmined/biovault/main/install.sh | bash
```

## TODO

- [x] CLI Tool
- [x] Data Format
- [ ] wizard
    - [x] create participant record
    - [x] checks for dependencies
    - [x] bv setup
    - [x] setup dependencies
        - [x] colab
        - [x] java
        - [x] nextflow
        - docker
    - [x] fetch mock data
    - submit analysis
    - [x] bv participant create
    - [x] bv participant list
    - [x] bv participant remove
- [x] bv project create
- [x] bv run ./project participant.yaml
- [x] bv biobank publish
- [x] bv biobank unpublish
- [x] bv biobank list
- [ ] bv update
- [ ] installer check existing install and version
- [ ] toggle docker mode
- [ ] include common modules like bcftools
- [ ] download deduplicate hashing and symlinking
- [ ] bv biobank list
  - show public path?
- [ ] windows symlinks plus tests

## Data Formats

- [ ] Change deep linking from:
  syft://madhava@openmined.org/private/biovault/participants.yaml#participants/MADHAVA
to:
  syft://madhava@openmined.org/private/biovault/participants.yaml#participants.MADHAVA


Participants in your biobank are kept in a private file like so:
`~/.biobank/participants.yaml`
```yaml
participants:
  MADHAVA:
    id: MADHAVA
    ref_version: GRCh38
    ref: /some/path/Homo_sapiens_assembly38.fasta
    ref_index: /some/path/Homo_sapiens_assembly38.fasta.fai
    aligned: /some/path/Madhava.cram
    aligned_index: /some/path/Madhava.cram.crai
```

When you choose to publish them an entry is added to a public syftbox file:
`~/SyftBox/datasites/madhava@openmined.org/public/biovault/participants.yaml`

```yaml
private_url: "syft://madhava@openmined.org/private/biovault/participants.yaml"

# Mock data anchors for testing
mock_data_grch38: &mock_data_grch38
  ref_version: GRCh38
  ref: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
  ref_index: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
  aligned: https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/data/ERR3239276/NA06985.final.cram
  aligned_index: https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/data/ERR3239276/NA06985.final.cram.crai
  # BLAKE3 checksums for verification (much faster than SHA-256)
  ref_b3sum: "49cbaceaf79ebc1da6581b2f7599cb03e6552ccce87584d1a0eaec59c3629368"
  ref_index_b3sum: "002cf8e0066a2226616b5d9cc09994ac06831cd907e13e521bef6dc69403d147"
  aligned_b3sum: "4556b84f32e58e1a5c4d7238352e9fc0bcaabd2478250733252f2b76047ba3ca"
  aligned_index_b3sum: "6914d3c6842670bdde272b8cc4dfaf858a84f379e9e79d8b24c1a89d577262e2"

participants:
  MADHAVA:
    id: MADHAVA
    url: "{root.private_url}#participants.MADHAVA"
    ref_version: GRCh38
    ref: "{url}.ref"
    ref_index: "{url}.ref_index"
    aligned: "{url}.aligned"
    aligned_index: "{url}.aligned_index"
    mock: *mock_data_grch38
```

This has the following schemas:

```yaml
---
schema: "org.openmined.biovault.participants-list-v1.0.0-beta.1"
title: "BioVault Participants List"
required: [participants]
properties:
  participants:
    type: map
    description: "Map of participant_id -> participant"
    key_types: string
    value_types: *participant

defs:
  participant: &participant
    type: object
    description: >
      Participant record. The `id` SHOULD be the same as the key used in the `participant` map.
    required: [id, ref_version, ref, ref_index, aligned, aligned_index]
    properties:
      id:
        type: string
        description: "Participant ID; should equal the key in the map"
      ref_version:
        type: enum
        values: ["GRCh38", "GRCh37"]
      ref:
        type: filepath
      ref_index:
        type: filepath
      aligned:
        type: filepath
      aligned_index:
        type: filepath

examples:
  - participants:
      TEST:
        id: TEST
        ref_version: GRCh38
        ref: ../data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
        ref_index: ../data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
        aligned: ../data/ERR3239283/NA07357.final.cram
        aligned_index: ../data/ERR3239283/NA07357.final.cram.crai

```

This can be added to your `resources.yaml`
```yaml
---
resources:
- name: biovault-participants
  path: syft://madhava@openmined.org/public/biovault/participants.yaml
  schema: org.openmined.biovault.participants-list-v1.0.0-beta.1
  schema_ref: ./resources/schemas/org.openmined.biovault.participants-list-v1.0.0-beta.1.yaml
```



## fastq combine
A utility to help combine large folders of fastq files.
- Checks for duplicates, mismatches and missing files
- Detects ONT metadata
- Suggests output filename
- Generates stats of input files using seqkit
- Verifies stats match on output file
- Creates blake3 hash file

Real Example:

```
bv fastq combine /Users/madhavajay/dev/carigenetics/onedrive/Fastqs/fastq_pass /Users/madhavajay/dev/carigenetics/combined
```

```
Scanning for FASTQ files in: /Users/madhavajay/dev/carigenetics/onedrive/Fastqs/fastq_pass

📁 File group: PBE09234_pass_2ce5fdd6_00ebc536
  Files: 433 (sequences 0-432)
  Total size: 23.69 GB

📁 File group: PBE09980_pass_290adfbe_f615b940
  Files: 378 (sequences 0-377)
  Total size: 83.84 GB

🔬 Checking for Oxford Nanopore metadata...

🧬 Oxford Nanopore Metadata Summary:
==================================================
  Run IDs found: 2
    - f615b940c5828ce0638002a4f0eb7f2ac23a6e16
    - 00ebc53629f8efdcca95471907a7c16cc0f4df12
  Flow cells: {"PBE09234", "PBE09980"}
  First start time: 2025-07-28 17:14:55.040747 -03:00
  Last start time: 2025-08-02 12:16:28.661769 -03:00
  Duration span: 115 hours 1 minutes
  Protocol groups: {"EntirelyYou"}
  Sample IDs: {"ENT0001"}
  Basecall models: {"dna_r10.4.1_e8.2_400bps_hac@v4.3.0"}

📊 Summary:
  Total files: 811
  Total size: 107.54 GB

❓ '/Users/madhavajay/dev/carigenetics/combined' doesn't exist and has no extension. Create as directory? (y/n): y
📁 Created output directory: /Users/madhavajay/dev/carigenetics/combined

📝 Suggested output filename: /Users/madhavajay/dev/carigenetics/combined/PBE09234-PBE09980-ENT0001-ONT-20250802-all.fastq.gz
Accept this filename? (y/n/e to edit): e
Enter new filename (without path) [PBE09234-PBE09980-ENT0001-ONT-20250802-all.fastq.gz]: MadhavaJay-WGS-carigenetics-PBE09234-PBE09980-ENT0001-ONT-20250802-all.fastq.gz
Using filename: /Users/madhavajay/dev/carigenetics/combined/MadhavaJay-WGS-carigenetics-PBE09234-PBE09980-ENT0001-ONT-20250802-all.fastq.gz

❓ Do you want to validate all files before combining? (y/n): y


```