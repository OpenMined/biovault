# Development Guide

## Prerequisites

- [just](https://github.com/casey/just) - Command runner
- Docker and Docker Compose
- Git with submodules support
- Rust and Cargo (for running integration tests)

## Setup

The project includes SyftBox as a submodule. After cloning:

```bash
git submodule update --init --recursive
```

## Available Commands

### Integration Testing

BioVault supports two integration testing modes:

#### Docker-based Testing (Containerized Clients)

Run the full integration test suite with Docker containers:
```bash
just test-integration-docker
```

Run tests and keep services running for inspection:
```bash
just test-integration-docker-inspect
# or
just test-integration-docker false
```

Clean up Docker test resources:
```bash
just cleanup-integration-test-docker
```

#### Local Testing (sbenv-based Clients)

Run the full integration test suite with local sbenv clients:
```bash
just test-integration-local
```

Run tests and keep services running for inspection:
```bash
just test-integration-local-inspect
# or
just test-integration-local false
```

Clean up local test resources:
```bash
just cleanup-integration-test-local
```

#### What the Tests Do

Both test modes:
- Check and configure `/etc/hosts` for MinIO
- Start SyftBox server with MinIO storage
- Start three SyftBox clients including a "bad actor" client for permission testing
- Run file synchronization tests between clients
- Test permission system with various access control scenarios
- Clean up automatically on success (unless disabled)
- Leave everything running on failure for debugging

The integrated permission tests verify:
- User-only read/write permissions
- Everyone read, user-only write permissions
- Specific user access controls
- Permission updates and propagation
- Bad actor access prevention

Test clients:
- **client1@syftbox.net**: Primary test client
- **client2@syftbox.net**: Secondary test client
- **bad@syftbox.net**: Bad actor client for security testing

#### Key Differences

- **Docker mode**: Uses `test-clients-docker/` directory, clients run in containers
- **Local mode**: Uses `test-clients-local/` directory, clients run as local processes via sbenv

### SyftBox Server Management

Start the SyftBox server with MinIO storage:
```bash
just start-syftbox-server
```

Stop the SyftBox server:
```bash
just stop-syftbox-server
```

Stop server and delete all volumes:
```bash
just stop-syftbox-server-delete-volumes
```

### SyftBox Client Management

#### Docker Clients

Start both Docker clients:
```bash
just start-syftbox-clients-docker
```

Start individual Docker clients:
```bash
just start-syftbox-client1-docker
just start-syftbox-client2-docker
```

Stop all Docker clients:
```bash
just stop-syftbox-clients-docker
```

#### Local Clients (sbenv)

Build sbenv tool:
```bash
just build-sbenv
```

Start both local clients:
```bash
just start-syftbox-clients-local
```

Start individual local clients:
```bash
just start-syftbox-client1-local
just start-syftbox-client2-local
```

Stop all local clients:
```bash
just stop-syftbox-clients-local
```

### Local Development Setup

Add MinIO hostname to /etc/hosts (required for clients running outside Docker):
```bash
just add-minio-host
```

This adds `127.0.0.1 minio` to your `/etc/hosts` file, allowing local clients to resolve the MinIO storage endpoint.

## Server Details

- **Server URL**: http://localhost:8080
- **MinIO Storage**: Accessible via the `minio` hostname when properly configured
- **Client 1**: http://localhost:7938 (client1@syftbox.net)
- **Client 2**: http://localhost:7939 (client2@syftbox.net)
- **Bad Client**: http://localhost:7940 (bad@syftbox.net) - Used for permission testing

## GitHub Workflows

The project runs nightly integration test workflows:

1. **Docker Integration Tests** (`biovault-integration-test-docker.yml`)
   - Runs at 2:00 AM UTC nightly
   - Tests with containerized clients
   - Includes permission system testing with bad actor client
   - Updates submodules to latest versions
   - Can be manually triggered with options

2. **Local Integration Tests** (`biovault-integration-test-local.yml`)
   - Runs at 2:30 AM UTC nightly
   - Tests with local sbenv clients
   - Includes permission system testing with bad actor client
   - Updates submodules to latest versions
   - Can be manually triggered with options

All workflows automatically update submodules to the latest versions during nightly runs to catch breaking changes in dependencies.

## Debugging

When tests fail or when using inspect mode, you can examine the running services:

### Docker Mode
```bash
# View server logs
docker logs syftbox-server

# View client logs
docker logs syftbox-client-client1-syftbox-net
docker logs syftbox-client-client2-syftbox-net

# Check synced files
ls -la test-clients-docker/

# Clean up when done
just cleanup-integration-test-docker
```

### Local Mode
```bash
# View server logs
docker logs syftbox-server

# View client logs
cat test-clients-local/client1@syftbox.net/.syftbox/logs/client.log
cat test-clients-local/client2@syftbox.net/.syftbox/logs/client.log

# Check synced files
ls -la test-clients-local/

# Clean up when done
just cleanup-integration-test-local
```