# BioVault

### Quick Install (One-liner)

```bash
curl -sSL https://raw.githubusercontent.com/openmined/biovault/main/install.sh | bash
```

- [x] CLI Tool
- [ ] data format
- [ ] wizard
    - create patient record
    - [x] checks for dependencies
    - [x] bv setup
    - [x] setup dependencies
        - [x] colab
        - [x] java
        - [x] nextflow
        - docker
    - [x] fetch mock data
    - submit analysis
    - bv patient create
    - bv patient list
    - bv patient remove

- [x] bv project create
- [x] bv run ./project patient.yaml
- [ ] bv update
- [ ] installer check existing install and version
- [ ] toggle docker mode
- [ ] include common modules like bcftools




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

```yaml
---
resources:
- name: biovault-participants
  path: syft://madhava@openmined.org/private/biovault/participants.yaml
  schema: org.openmined.biovault.participants-list-v1.0.0-beta.1
  schema_ref: ./resources/schemas/org.openmined.biovault.participants-list-v1.0.0-beta.1.yaml


participants:
  TEST:
    id: TEST
    ref_version: GRCh38
    ref: ../data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
    ref_index: ../data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
    aligned: ../data/ERR3239283/NA07357.final.cram
    aligned_index: ../data/ERR3239283/NA07357.final.cram.crai
```

private_url: `syft://madhava@openmined.org/private/biovault/participants.yaml#/participants/TEST`

syft://{datasite}/public/biovault/participants.yaml

``yaml
---
schema: "net.syftbox.twin-file-v1.0.0-beta.1"
name: "participants"
description_path: "syft://{datasite}/public/biovault/participants/DATA.md"
format: "yaml"
real_path: syft://{datasite}/private/biovault/participants.yaml
mock_path: syft://{datasite}/public/biovault/participants.mock.yaml
```
