# BioVault Hello World Tutorial

Welcome to your first BioVault tutorial! This guide will walk you through setting up BioVault and running your first genomic analysis.

## What is BioVault?

BioVault is a free, open-source platform for collaborative genomics that allows researchers to share insights without ever sharing raw genetic data. It uses end-to-end encryption and secure enclaves to protect privacy while enabling scientific collaboration.

## What You'll Learn

In this tutorial, you will:
- Install BioVault and its dependencies
- Set up a complete SNP counting project with 23andMe data
- Run a simple analysis to count genetic variants
- Process over 600,000 SNPs from a real genome file
- Understand the basic structure of a BioVault project

## Prerequisites

Before starting, make sure you have:
- A computer with internet access
- Basic familiarity with command line/terminal
- About 30 minutes to complete the tutorial
- At least 100MB of free disk space

**Note:** This tutorial uses Y chromosome data because it's the smallest chromosome, making it perfect for learning without requiring massive downloads. We'll build a SNP counting analysis that can work with various genetic data formats.

## Step 1: Install BioVault

The BioVault command line tool is called `bv`. You have several installation options:

### Option A: Quick Install (Recommended)
This is the easiest way to install BioVault:
```bash
curl -sSL https://raw.githubusercontent.com/openmined/biovault/main/install.sh | bash
```

### Option B: Manual Download
If you prefer to download manually, visit: https://github.com/openmined/biovault

### Option C: Using Cargo (for Rust developers)
If you have Rust installed:
```bash
cargo install biovault
```

## Step 2: Update BioVault
Make sure you have the latest version:
```bash
bv update
```


## Step 3: Check Dependencies
BioVault requires several dependencies to work properly. Let's check if your system has everything needed:
```bash
bv check
```

You should see output like this:
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

### Important: First-time SyftBox Setup

When you install SyftBox for the first time, you'll need to create a special ignore file to ensure proper syncing:

1. Navigate to the datasites directory:
   ```bash
   cd ~/Syftbox/datasites
   ```

2. Create a file called `syftignore` with the following contents:
   ```
   **/app_data/slack-mcp/**
   khoa@openmined.org
   koen@openmined.org
   amita.j.shukla@gmail.com
   dhingra.atul92@gmail.com
   zach@empire.email
   ```

3. Restart SyftBox:
   ```bash
   syftbox
   ```

**Note:** Keep SyftBox running in this terminal and open a new terminal for additional work.

### If You Have Missing Dependencies
If any dependencies are missing, you can try the automatic setup:
```bash
bv setup
```

This command will attempt to install missing dependencies on supported systems (macOS and Google Colab).

## Step 4: Initialize BioVault
Set up BioVault for your user account:
```bash
bv init
```

This command will ask for your email address. **Important:** Make sure this matches your SyftBox email if you have one.

**Note:** The first thing you want to do is test your code on sample data and run it locally to make sure your code works before sending your analysis to someone else.
## Step 5: Set Up Sample Data
For this tutorial, we'll work with a 23andMe genome file format. This is a simple text file containing genetic variants (SNPs) that's perfect for learning.

Let's copy the count-snps example which includes a sample 23andMe genome file:
```bash
cp -r /Users/dawnxi/biovault/cli/examples/count-snps ./count-snps
```

This copies:
- A complete count-snps project template
- A sample 23andMe genome file (`genome_23andMe_v4_Full.txt`) with over 600,000 genetic variants
- All necessary workflow and configuration files

**Note:** The 23andMe file format is a simple text file with genetic variants, making it perfect for learning SNP counting without requiring massive downloads or complex genomic data processing.

## Step 6: Explore the Project Structure
Let's look at what we copied and understand the project structure:
```bash
cd ./count-snps
ls -la
```

You should see:
```
count-snps/
├── assets/
│   └── count_number_of_snps.py
├── genome_23andMe_v4_Full.txt
├── project.yaml
├── template.nf
└── workflow.nf
```

This is a complete BioVault project ready for SNP counting analysis!

## Step 7: Understand the Project Structure
Let's examine the key files in our count-snps project:

### The 23andMe Genome File
```bash
head -5 genome_23andMe_v4_Full.txt
```

You'll see the format:
```
# rsid	chromosome	position	genotype
rs12564807	1	734462	AA
rs3131972	1	752721	AG
rs148828841	1	760998	AC
rs12124819	1	776546	AA
```

This file contains over 600,000 genetic variants (SNPs) from a 23andMe genome test.

### Project Files
- **`project.yaml`**: Contains project settings, participant information, and analysis parameters
- **`workflow.nf`**: The main analysis workflow that processes the 23andMe data
- **`template.nf`**: The Nextflow template that sets up the workflow execution
- **`assets/count_number_of_snps.py`**: Python script that counts SNPs in the input data

### How It Works
The workflow will:
1. Read the 23andMe genome file
2. Count the number of genetic variants (SNPs)
3. Output the total count to a results file

## Step 8: Run Your First Analysis
Now let's run the SNP counting analysis on our 23andMe genome file:
```bash
bv run ./count-snps genome_23andMe_v4_Full.txt
```

This command tells BioVault to:
1. Use the project in the `./count-snps` directory
2. Analyze the 23andMe genome file `genome_23andMe_v4_Full.txt`
3. Count the number of SNPs (genetic variants) in the data

### What to Expect
You should see output similar to this:
```
Launching `~/.biovault/env/default/template.nf`

executor >  local (1)
Participant ID: genome_23andMe_v4_Full
SNP file: /path/to/genome_23andMe_v4_Full.txt
Assets Directory: /path/to/assets
Results Directory: /path/to/results

[abc12345] USER:count_number_of_snps (1) | 1 of 1 ✔

===== Number of SNPs: 601804 =====

Workflow completed successfully!
```

The analysis will:
1. Process the 23andMe genome file
2. Count the number of genetic variants (SNPs)
3. Display the total count (should be around 601,804 SNPs)
4. Save the results to a file

## Step 9: View Your Results
Check the results of your analysis:
```bash
ls ./count-snps/results
cat ./count-snps/results/number_of_snps.txt
```

This will show you the output files generated by your SNP counting analysis, including the total number of SNPs found in the 23andMe genome file.

## 🎉 Congratulations!
You've successfully:
- Installed and configured BioVault
- Set up a complete SNP counting project
- Run your first genomic analysis on 23andMe data
- Counted over 600,000 genetic variants (SNPs)
- Learned how BioVault processes genetic data files

## What's Next?
- **Learn more**: Check out [Tutorial 2: Submit Your Project](tutorials/2_submit_your_project.md) to learn how to share your analysis with others
- **Create a biobank**: See [Tutorial 3: Create a Biobank](tutorials/3_create_biobank.md) to set up your own data sharing network
- **Explore**: Try modifying the `workflow.nf` file to customize your analysis or work with different types of genetic data

## Troubleshooting
If you encounter any issues:
1. Make sure all dependencies are installed: `bv check`
2. Try the automatic setup: `bv setup`
3. Check the [BioVault documentation](https://biovault.net) for more help
