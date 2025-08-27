# BioVault Google Colab Test Environment

This directory contains Docker-based tests to emulate the Google Colab environment for testing BioVault installation and setup procedures.

## Environment Details

The test container emulates:
- **OS**: Ubuntu 22.04.4 LTS (jammy)
- **Architecture**: x86_64
- **Java**: OpenJDK 11.0.28 (default Colab version)
- **Environment**: `COLAB_RELEASE_TAG=release-colab_20250822-060052_RC00`
- **Build Tools**: Rust toolchain for building `bv` binary

## Quick Start

### Run All Tests
```bash
./test_colab_install.sh
```

This script will:
1. Build the Docker container
2. Mount the source code
3. Build the `bv` binary inside the container
4. Run `bv setup` to install dependencies
5. Verify installations

### Manual Testing with Docker Compose

```bash
# Build and start the container
docker-compose up -d --build

# Enter the container
docker-compose exec colab-test bash

# Inside the container:
cp -r /workspace/cli_src /workspace/cli
cd /workspace/cli
cargo build --release
./target/release/bv setup
./target/release/bv check

# Clean up
docker-compose down
```

### Manual Testing with Docker

```bash
# Build the image
docker build -t biovault-colab-test .

# Run interactively
docker run --rm -it \
  -v "$(pwd)/../../cli:/workspace/cli:ro" \
  biovault-colab-test

# Inside the container:
cp -r /workspace/cli_src /workspace/cli
cd /workspace/cli
cargo build --release
./target/release/bv setup
./target/release/bv check
```

## What Gets Tested

The setup process will:
1. **Detect** the Google Colab environment via `COLAB_RELEASE_TAG`
2. **Upgrade Java** from 11 to 21 using `apt-get`
3. **Skip Docker** (not supported in Colab)
4. **Install Nextflow** to `/usr/local/bin`
5. **Verify** all installations

## Expected Results

After running `bv setup`:
- Java should be upgraded to version 21
- Nextflow should be installed and accessible
- Docker will be skipped with an explanation
- `bv check` should show all dependencies satisfied (except Docker)