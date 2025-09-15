# BioVault Hello World

- In this tutorial we will:
- install biovault and its dependencies
- download a small 200mb slice of someones genome with just the Y Chromosome (because its the smallest chromosome)
- run a haplogroup test against this chromosome to do some really basic ancestry estimation

1) install bv
Curl:
```
curl -sSL https://raw.githubusercontent.com/openmined/biovault/main/install.sh | bash
```

Manual Download:
https://github.com/openmined/biovault

Cargo:
```
cargo install biovault
```

2) make sure bv is up to date
```
bv update
```

3) check your system has the correct depenencies
```
bv check
```

```
BioVault Dependency Check
=========================

Checking java...  (version 23)✓ Found
Checking docker... ✓ Found (running)
Checking nextflow... ✓ Found
Checking syftbox... ✓ Found

=========================
✓ All dependencies satisfied!
```

4) Missing Dependencies
If you have any missing dependencies follow the instructions or run:
```
bv setup
```

5) Initialize
Run `bv init` which will ask for your email.
Make sure this matches your syftbox email.

6) Sample Data
Lets get some sample data on our machine, we will download a Y Chromosome which is only 100mb and its reference which is 50mb.
`bv sample-data fetch NA06985-chrY`.

This data will be stored here: `~/.biovault/data/sample`

7) Create example project
```
bv project create --haplo-y
```

8) open in your favourite editor
code ./haplo-y

9) The format of a project is as follows:
```
.
├── assets <- optional assets
├── project.yaml <- required
└── workflow.nf <- default entry point
```

10) run the example on your sample
```
bv run ./haplo-y NA06985-chrY
```

You should see:
```
Launching `~/.biovault/env/default/template.nf`

executor >  local (5)
[3f/91db27] USE…lect_panel (select_panel) | 1 of 1 ✔
[68/c9b331] USER:build_regions (regions)  | 1 of 1 ✔
[41/198880] USER:call_sites (1)           | 1 of 1 ✔
[7e/90870d] USER:query_to_table (1)       | 1 of 1 ✔
[38/d8760f] USER:interpret_haplogroup (1) | 1 of 1 ✔

Workflow completed successfully!
```

The results will be in:
ls ./haplo-y/results

Congratulations you just created and ran your first genomic analysis of the Y Chromosome.

Next: How to submit your project to another biovault biobank - [2_submit_your_project.md](tutorials/2_submit_your_project.md).