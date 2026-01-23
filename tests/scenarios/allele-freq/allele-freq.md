# syqure-allele-freq scenario plan

## Goal
Create a new scenario called `syqure-allele-freq` that:
- Brings up 3 datasites (aggregator, client1, client2).
- Imports and submits the allele-freq project from the aggregator to both clients.
- Generates 100 synthetic genotype files per client in parallel:
  - client1 with elevated APOL1 variants.
  - client2 with elevated THALASSEMIA variants.
- Runs allele frequency computation first, producing `allele_freq.npz` + `locus_index.json` per datasite.
- Runs a Syqure step that securely sums the allele-frequency arrays across parties and returns the final sum.

## Existing pieces to reuse
- `biovault/tests/scenarios/allele-freq/allele-freq/` (project + workflow)
  - Outputs: `allele_freq.tsv`, `allele_freq.npz`, `locus_index.json`.
- `biovault/tests/scenarios/allele-freq/test-allele-freq.sh`
  - Shows how to generate synthetic genotype files with APOL1/THAL overlays.
- `biovault/tests/scenarios/syqure-project.yaml`
  - Pattern for submit + inbox + process of a Syqure project across 3 datasites.
- `biovault/tests/scenarios/syqure/aggregate.sh` + `syqure-aggregate.yaml`
  - File-transport Syqure helper flow (ACLs, waiting, local/datasite paths).

## Proposed new files
1) Scenario runner:
   - `biovault/tests/scenarios/syqure-allele-freq.yaml`

2) Project directory:
   - `biovault/tests/scenarios/allele-freq/syqure-allele-freq/`
     - `project.yaml`
     - `workflow.nf` (or a shell workflow that wraps Nextflow + Syqure)
     - `pipeline.yaml` (optional if reusing `bv run` patterns)
     - `assets/` (Syqure codon program + optional template inputs)

3) Data generation helper (new script):
   - `biovault/tests/scripts/gen_allele_freq_data.sh`
     - Params: output_dir, count, seed, overlay variants JSON.
     - Uses `bvs synthetic` to emit 100 genotype files.
     - Produces a samplesheet CSV (participant_id, genotype_file).

## Scenario flow (syqure-allele-freq.yaml)
High-level step list:

1) Devstack setup
- Start 3 clients with `devstack.sh --reset`.

2) Ensure datasite dirs exist
- Create datasites root under each sandbox.

3) Submit allele-freq project
- Copy project directory to aggregator datasite.
- `bv submit` from aggregator to aggregator/client1/client2 (same as syqure-project).

4) Generate genotype data in parallel
- On client1: generate 100 genotype files with APOL1 overlay at high frequency.
- On client2: generate 100 genotype files with THAL overlay at high frequency.
- Optionally run in parallel using background jobs or separate scenario steps.

Suggested frequencies:
- APOL1: 0.6-0.8
- THAL: 0.6-0.8

5) Catalog or samplesheet creation per datasite
Two options (pick one):
- **Option A**: Build samplesheet directly from the output dir.
  - `bv samplesheet create <data_dir> <samplesheet.csv> --file_filter "*.txt"`
- **Option B**: Catalog files first, then export samplesheet.
  - `bv files import <data_dir> --pattern "{participant_id}.txt" --non-interactive`
  - `bv samplesheet export-catalog <samplesheet.csv> --data-type genotype`

6) Process the project messages
- Same pattern as `syqure-project.yaml`:
  - sync inbox on all 3 datasites
  - extract message IDs
  - run `tests/scripts/run_syqure_project.sh` (or a new runner if needed)

7) Workflow behavior inside the project
- Step A: allele-freq computation
  - For each datasite, run allele-freq workflow using its samplesheet.
  - Outputs `allele_freq.npz` + `locus_index.json` per datasite.

- Step B: align indices for secure sum
  - Need a shared, consistent locus index for all parties.
  - Use `allele_freq_utils.py` to build a union index and realign each datasite’s arrays.
  - Write aligned `*.npz` for Syqure input.

- Step C: Syqure secure sum
  - Use a Syqure codon program that loads the aligned arrays and computes element-wise sum.
  - Output final summed array as a file in results.

## Notes on the Syqure codon step
- Codon must read numpy arrays (`.npz`) or a simpler binary/TSV encoding.
- If Codon cannot read `.npz` directly, add a lightweight conversion step to CSV/flat binary before Syqure.
- Output could be:
  - a summed `.npz`
  - a TSV with rsid + total AC/AN
  - an aggregate summary for APOL1/THAL loci as a quick sanity check.

## Open questions / decisions needed
1) **Input wiring**: how do we pass each datasite’s samplesheet to the project?
   - If `bv message process` supports input prompts, we may need to prepopulate or pass via env.
2) **Union index generation**: should aggregator compute and distribute a union `locus_index.json`?
3) **Syqure I/O format**: does the Syqure Codon environment support reading `.npz`?
4) **Where to run Syqure**: inside Nextflow (workflow.nf) or in a shell workflow wrapper?

## Next steps to implement (once decisions above are confirmed)
- Create `syqure-allele-freq.yaml` scenario.
- Add `gen_allele_freq_data.sh` helper for deterministic data generation + samplesheet creation.
- Add `syqure-allele-freq` project folder with:
  - workflow updated to accept samplesheet input
  - optional index alignment helper
  - Syqure aggregation step + codon program
- Wire scenario to submit + process project and verify results for APOL1/THAL rsids.
