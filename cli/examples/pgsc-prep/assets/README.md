# PGS Catalog Preparation (pgsc-prep)

This workflow converts raw SNP genotype files (like `000000_carika.txt`) into plink binary format, which is required for calculating polygenic scores using `pgsc_calc`.

## What This Workflow Does

1. **Reads SNP file** - Parses tab-separated genotype data
2. **Converts to plink format** - Creates `.ped` and `.map` text files
3. **Creates binary files** - Generates `.bed`, `.bim`, `.fam` files (more efficient)
4. **Outputs ready-to-use files** - All files needed for pgsc_calc

## Output Files

After running this workflow, you'll have:

- `{participant_id}.bed` - Binary genotype file
- `{participant_id}.bim` - Variant information (SNP IDs, positions, alleles)
- `{participant_id}.fam` - Sample information

These files are in plink binary format and can be used directly with pgsc_calc.

## How to Run pgsc_calc

Once you have the plink files, you can calculate polygenic scores using the PGS Catalog calculator:

### Step 1: Install/Setup pgsc_calc

You'll need Nextflow and Docker installed. Then:

```bash
# No installation needed - nextflow will fetch it automatically
# Just make sure you have Docker running
```

### Step 2: Run pgsc_calc

Basic example to calculate a single PGS:

```bash
nextflow run pgscatalog/pgsc_calc \
  -profile docker \
  --input /path/to/participant.bed \
  --pgs_id PGS000001 \
  --outdir ./pgsc_results
```

### Step 3: Calculate Multiple PGS

You can calculate multiple polygenic scores at once:

```bash
nextflow run pgscatalog/pgsc_calc \
  -profile docker \
  --input /path/to/participant.bed \
  --pgs_id PGS000001,PGS000002,PGS000003 \
  --outdir ./pgsc_results
```

### Step 4: Use a Scorefile

If you have a custom scoring file:

```bash
nextflow run pgscatalog/pgsc_calc \
  -profile docker \
  --input /path/to/participant.bed \
  --scorefile /path/to/custom_score.txt \
  --outdir ./pgsc_results
```

## Common PGS Catalog Scores

Here are some popular polygenic scores you might want to calculate:

- **PGS000001** - Type 2 Diabetes
- **PGS000002** - Coronary Artery Disease
- **PGS000018** - Breast Cancer
- **PGS000027** - Height
- **PGS000039** - Body Mass Index (BMI)
- **PGS000662** - Alzheimer's Disease

Browse more scores at: https://www.pgscatalog.org/

## Advanced Options

### Ancestry Adjustment

If you want to adjust scores based on genetic ancestry:

```bash
nextflow run pgscatalog/pgsc_calc \
  -profile docker \
  --input /path/to/participant.bed \
  --pgs_id PGS000001 \
  --run_ancestry reference_panel \
  --outdir ./pgsc_results
```

### Specify Genome Build

```bash
nextflow run pgscatalog/pgsc_calc \
  -profile docker \
  --input /path/to/participant.bed \
  --pgs_id PGS000001 \
  --target_build GRCh38 \
  --outdir ./pgsc_results
```

## Output from pgsc_calc

After pgsc_calc completes, you'll find:

- `aggregated_scores.txt` - Final polygenic scores
- `score_distributions.html` - Interactive visualization
- `matches/` - Details about which SNPs matched
- `logs/` - Execution logs

## Troubleshooting

### Low Match Rate

If you get a low match rate warning:

- Your SNP file might use a different genome build (GRCh37 vs GRCh38)
- Some SNPs might have different IDs (rs numbers changed between builds)
- Solution: Specify `--target_build` parameter

### Memory Issues

If you run out of memory:

- Add `-profile docker,test` to use smaller test resources
- Increase Docker memory limits in Docker Desktop settings

### Docker Issues

If Docker fails:

- Make sure Docker Desktop is running
- Try `-profile singularity` instead if available

## More Information

- **pgsc_calc documentation**: https://pgsc-calc.readthedocs.io/
- **PGS Catalog**: https://www.pgscatalog.org/
- **BioVault docs**: https://github.com/OpenMined/biovault

## Example Complete Workflow

```bash
# 1. Run this BioVault workflow to convert format
bv run --project pgsc-prep --participant 000000

# 2. Navigate to results directory
cd results-real/000000/

# 3. Run pgsc_calc
nextflow run pgscatalog/pgsc_calc \
  -profile docker \
  --input 000000.bed \
  --pgs_id PGS000001,PGS000002 \
  --outdir ./pgs_results

# 4. View results
cat pgs_results/aggregated_scores.txt
open pgs_results/score_distributions.html
```

## Notes

- This workflow assumes your SNP file is in GRCh38 format (most modern files are)
- If your data is GRCh37, specify `--target_build GRCh37` when running pgsc_calc
- The conversion process preserves all SNPs but pgsc_calc will only use those in the scoring file
- Large SNP files (600k+ variants) may take 5-10 minutes to process
