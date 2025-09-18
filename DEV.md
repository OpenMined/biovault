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

Run the full integration test suite (default - cleans up on success):
```bash
just test-integration-docker
```

Run tests and keep services running for inspection:
```bash
just test-integration-docker-inspect
# or
just test-integration-docker false
```

Clean up test resources after inspection:
```bash
just cleanup-integration-test
```

The integration tests:
- Check and configure `/etc/hosts` for MinIO
- Start SyftBox server with MinIO storage
- Start two SyftBox client containers
- Run file synchronization tests between clients
- Clean up automatically on success (unless disabled)
- Leave everything running on failure for debugging

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

Start both test clients:
```bash
just start-syftbox-clients
```

Start individual clients:
```bash
just start-syftbox-client1
just start-syftbox-client2
```

Stop all clients:
```bash
just stop-syftbox-clients
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

## Debugging

When tests fail or when using `test-integration-docker-inspect`, you can inspect the running services:

```bash
# View server logs
docker logs syftbox-server

# View client logs
docker logs syftbox-client-client1-syftbox-net
docker logs syftbox-client-client2-syftbox-net

# Check synced files
ls -la test-clients/

# Clean up when done
just cleanup-integration-test
```