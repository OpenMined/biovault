# Pipeline & SQL Store Enhancements

- **Top-Level Inputs**: Pipelines now declare reusable inputs once (e.g. `samplesheet`, `data_dir`), and steps reference them via `inputs.<name>`. Validation lists these inputs first, and `bv run` honors overrides such as `--set inputs.samplesheet=/path/to.csv`.

- **Sequential Runner Updates**: `bv run pipeline.yaml` executes every step in order, handling pipeline/step overrides, `--results-dir`, and automatically wiring outputs for downstream steps. (The sandbox fails due to missing Java, but on a machine with a JRE the full run completes.)

- **SQL Store Support**: A step can add `store: …` with `kind: sql`, `destination: SQL()` (defaults to the BioVault SQLite DB), `source`, `table_name`, and `key_column`. After a step finishes, its CSV/TSV is ingested into the DB; the CLI logs the source file, destination table, and target database path.

- **Example Pipeline**: `pipeline_sql.yaml` chains the sample projects (`filter-samples → count-lines → sum-total`) and writes the counted-sheet output to `pipeline_counts_{run_id}`. Overriding both inputs and running produces:

  ```bash
  ./bv run pipeline_sql.yaml \
    --set inputs.samplesheet=/…/pipeline/data/participants.csv \
    --set inputs.data_dir=/…/pipeline/data \
    --results-dir /…/results/demo-sql
  ```

  The `counts_sql` store logs a line such as:

  ```text
  💾 Stored 'counts_sql' output 'counted_sheet' into table pipeline_counts_20251023050324 (rows: …).
      source: /…/results/demo-sql/count/line_counts.csv
      database: /…/Desktop/BioVault/biovault.db
  ```

- **Workflow Script**: `run_demo_sql.sh` demonstrates the end-to-end flow—validate and then run with hard-coded paths to the sample data.

- **Documentation & Wizard**: Both pipeline docs reflect the new `inputs:` block and SQL store syntax; the wizard emits `inputs.*` bindings by default.

**Next steps**: Install Java to run Nextflow locally; once the environment has a JRE, re-run `./run_demo_sql.sh` to see the store in action and query the resulting table with `bv sql`.
