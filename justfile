DOCKER_DIR := "./syftbox/docker"
SYFTBOX_DIR := "./syftbox"
SBENV_DIR := "./sbenv"
TEST_CLIENTS_DOCKER_DIR := "./test-clients-docker"
TEST_CLIENTS_LOCAL_DIR := "./test-clients-local"

# Test client configurations
TEST_CLIENT1_EMAIL := "client1@syftbox.net"
TEST_CLIENT2_EMAIL := "client2@syftbox.net"
TEST_BAD_CLIENT_EMAIL := "bad@syftbox.net"
TEST_CLIENT1_NAME := "client1-syftbox-net"
TEST_CLIENT2_NAME := "client2-syftbox-net"
TEST_BAD_CLIENT_NAME := "bad-syftbox-net"
TEST_CLIENT1_PORT := "7938"
TEST_CLIENT2_PORT := "7939"
TEST_BAD_CLIENT_PORT := "7940"
TEST_SERVER_URL := "http://localhost:8080"

# Check if server is running
server-is-running:
    @docker ps --format 'table {{{{.Names}}}}' | grep -q syftbox-server && echo "true" || echo "false"

# Check if MinIO is running
minio-is-running:
    @docker ps --format 'table {{{{.Names}}}}' | grep -q syftbox-minio && echo "true" || echo "false"

# Reset MinIO data (keeps container running)
reset-minio:
    @echo "Resetting MinIO data..."
    @docker exec syftbox-minio sh -c "rm -rf /data/*" 2>/dev/null || true
    @echo "MinIO data cleared"

# Reset server (restart with clean state)
reset-server:
    @echo "Resetting SyftBox server..."
    @docker restart syftbox-server 2>/dev/null || true
    @echo "Waiting for server to restart..."
    @sleep 3
    @echo "Server reset complete"

start-syftbox-server:
    @if [ "$(just server-is-running)" = "true" ]; then \
        echo "SyftBox server already running"; \
    else \
        echo "Starting SyftBox server with MinIO..."; \
        cd {{DOCKER_DIR}} && COMPOSE_BAKE=true docker compose up -d --build minio server; \
        echo "Waiting for server to be ready..."; \
        sleep 5; \
        echo "Server started at http://localhost:8080"; \
    fi

stop-syftbox-server:
    cd {{DOCKER_DIR}} && docker compose down

stop-syftbox-server-delete-volumes:
    cd {{DOCKER_DIR}} && docker compose down -v

add-minio-host:
    @echo "Adding minio to /etc/hosts..."
    @if ! grep -q "^[^#]*127.0.0.1.*minio" /etc/hosts; then \
        echo "127.0.0.1 minio" | sudo tee -a /etc/hosts > /dev/null; \
        echo "Added minio to /etc/hosts"; \
    else \
        echo "minio already in /etc/hosts (not commented)"; \
    fi

# Generic docker client start command
start-syftbox-client-docker email name port:
    @echo "Starting docker client {{email}}..."
    @mkdir -p {{TEST_CLIENTS_DOCKER_DIR}}
    @echo "Building client image..."
    @if [ -n "${DOCKER_BUILDX:-}" ]; then \
        cd {{SYFTBOX_DIR}} && docker buildx build --cache-from=type=gha --cache-to=type=gha,mode=max -f docker/Dockerfile.client -t syftbox-client --load .; \
    else \
        cd {{SYFTBOX_DIR}} && docker build -f docker/Dockerfile.client -t syftbox-client .; \
    fi
    @echo "Starting client container..."
    docker run -d \
        --name syftbox-client-{{name}} \
        --network docker_syftbox-network \
        -p {{port}}:7938 \
        -e SYFTBOX_SERVER_URL=http://syftbox-server:8080 \
        -e SYFTBOX_AUTH_ENABLED=0 \
        -v "$(pwd)/{{TEST_CLIENTS_DOCKER_DIR}}:/data/clients" \
        syftbox-client {{email}}
    @echo "Client started at http://localhost:{{port}}"

# Start docker client 1
start-syftbox-client1-docker:
    just start-syftbox-client-docker "{{TEST_CLIENT1_EMAIL}}" "{{TEST_CLIENT1_NAME}}" "{{TEST_CLIENT1_PORT}}"

# Start docker client 2
start-syftbox-client2-docker:
    just start-syftbox-client-docker "{{TEST_CLIENT2_EMAIL}}" "{{TEST_CLIENT2_NAME}}" "{{TEST_CLIENT2_PORT}}"

# Start docker bad client
start-syftbox-bad-client-docker:
    just start-syftbox-client-docker "{{TEST_BAD_CLIENT_EMAIL}}" "{{TEST_BAD_CLIENT_NAME}}" "{{TEST_BAD_CLIENT_PORT}}"

# Start all docker clients (including bad client)
start-syftbox-all-clients-docker:
    @echo "Cleaning up old test data..."
    @if [ -n "${CI:-}" ] && [ -d {{TEST_CLIENTS_DOCKER_DIR}} ]; then \
        echo "In CI environment, using sudo to remove root-owned files"; \
        sudo rm -rf {{TEST_CLIENTS_DOCKER_DIR}}; \
    else \
        rm -rf {{TEST_CLIENTS_DOCKER_DIR}}; \
    fi
    @just start-syftbox-client1-docker
    @sleep 2
    @just start-syftbox-client2-docker
    @sleep 2
    @just start-syftbox-bad-client-docker

# Stop all docker clients
stop-syftbox-all-clients-docker:
    -docker stop syftbox-client-{{TEST_CLIENT1_NAME}}
    -docker stop syftbox-client-{{TEST_CLIENT2_NAME}}
    -docker stop syftbox-client-{{TEST_BAD_CLIENT_NAME}}
    -docker rm syftbox-client-{{TEST_CLIENT1_NAME}}
    -docker rm syftbox-client-{{TEST_CLIENT2_NAME}}
    -docker rm syftbox-client-{{TEST_BAD_CLIENT_NAME}}

# Wait for server to be ready
wait-for-server:
    @echo "Waiting for server to be ready..."
    @timeout 60 bash -c 'until curl -f http://localhost:8080 2>/dev/null; do echo "Waiting..."; sleep 2; done'
    @echo "Server is ready!"

# Quick reset for reruns - keeps docker stack, just cleans data
quick-reset:
    @echo "=== Quick Reset for Test Reruns ==="
    @if [ "$(just server-is-running)" = "false" ] || [ "$(just minio-is-running)" = "false" ]; then \
        echo "Server/MinIO not running, starting fresh..."; \
        just start-syftbox-server; \
    else \
        echo "Server/MinIO already running, resetting..."; \
        just reset-minio; \
        just reset-server; \
    fi
    @echo "Cleaning test client directories..."
    @echo "  Stopping any running local clients first..."
    @-just stop-syftbox-all-clients-local 2>/dev/null || true
    @echo "  Removing test directories..."
    @rm -rf {{TEST_CLIENTS_DOCKER_DIR}} 2>/dev/null || true
    @rm -rf {{TEST_CLIENTS_LOCAL_DIR}} 2>/dev/null || true
    @echo "Quick reset complete!"


# Run all integration tests with docker clients
run-integration-tests-docker:
    @echo "Building BioVault binary..."
    cd cli && cargo build --release
    @echo "Running integration tests with docker clients..."
    cd cli && \
    SYFTBOX_SERVER_URL=http://localhost:8080 \
    SYFTBOX_CLIENT1_EMAIL={{TEST_CLIENT1_EMAIL}} \
    SYFTBOX_CLIENT2_EMAIL={{TEST_CLIENT2_EMAIL}} \
    SYFTBOX_BAD_CLIENT_EMAIL={{TEST_BAD_CLIENT_EMAIL}} \
    TEST_CLIENTS_DIR=../{{TEST_CLIENTS_DOCKER_DIR}} \
    TEST_MODE=docker \
    cargo test --features e2e-tests --test '*' -- --ignored --nocapture

# Full integration test with Docker - starts everything, runs tests, cleans up on success
test-integration-docker cleanup="true" reset="false":
    @echo "=== Starting Integration Test with Docker ==="
    @if [ "{{reset}}" = "true" ] || [ "{{reset}}" = "--reset" ]; then \
        echo "Reset requested, cleaning data..."; \
        just quick-reset; \
    fi
    @echo "1. Checking /etc/hosts for minio..."
    @if ! grep -q "^[^#]*127.0.0.1.*minio" /etc/hosts; then \
        echo "❌ minio not found in /etc/hosts. Running: just add-minio-host"; \
        just add-minio-host; \
    else \
        echo "✓ minio found in /etc/hosts"; \
    fi
    @echo ""
    @echo "2. Starting SyftBox server..."
    @just start-syftbox-server
    @echo ""
    @echo "3. Waiting for server..."
    @just wait-for-server
    @echo ""
    @echo "4. Starting SyftBox docker clients (including bad client)..."
    @just start-syftbox-all-clients-docker
    @echo ""
    @echo "5. Waiting for clients to initialize..."
    @sleep 10
    @echo ""
    @echo "6. Running integration tests..."
    @if just run-integration-tests-docker; then \
        echo ""; \
        if [ "{{cleanup}}" = "true" ]; then \
            echo "✅ All tests passed! Cleaning up..."; \
            just stop-syftbox-all-clients-docker; \
            just stop-syftbox-server; \
            rm -rf {{TEST_CLIENTS_DOCKER_DIR}}; \
            echo "✅ Cleanup complete!"; \
        else \
            echo "✅ All tests passed! Leaving services running for inspection."; \
            echo ""; \
            echo "Services are still running. You can:"; \
            echo "  - Check server logs: docker logs syftbox-server"; \
            echo "  - Check client1 logs: docker logs syftbox-client-{{TEST_CLIENT1_NAME}}"; \
            echo "  - Check client2 logs: docker logs syftbox-client-{{TEST_CLIENT2_NAME}}"; \
            echo "  - Check bad client logs: docker logs syftbox-client-{{TEST_BAD_CLIENT_NAME}}"; \
            echo "  - Explore files: ls -la {{TEST_CLIENTS_DOCKER_DIR}}/"; \
            echo ""; \
            echo "When done, run:"; \
            echo "  just cleanup-integration-test-docker"; \
        fi \
    else \
        echo ""; \
        echo "❌ Tests failed! Leaving services running for debugging."; \
        echo ""; \
        echo "Debug commands:"; \
        echo "  - Server logs: docker logs syftbox-server"; \
        echo "  - Client1 logs: docker logs syftbox-client-{{TEST_CLIENT1_NAME}}"; \
        echo "  - Client2 logs: docker logs syftbox-client-{{TEST_CLIENT2_NAME}}"; \
        echo "  - Bad client logs: docker logs syftbox-client-{{TEST_BAD_CLIENT_NAME}}"; \
        echo "  - Check files: ls -la {{TEST_CLIENTS_DOCKER_DIR}}/"; \
        echo ""; \
        echo "When done debugging, run:"; \
        echo "  just cleanup-integration-test"; \
        exit 1; \
    fi

# Run integration test without cleanup (keeps services running after)
test-integration-docker-inspect:
    just test-integration-docker false false

# Run integration test with reset (clears data but keeps containers)
test-integration-docker-reset:
    just test-integration-docker true --reset


# Clean up docker integration test resources
cleanup-integration-test-docker:
    @echo "Cleaning up docker integration test resources..."
    -just stop-syftbox-all-clients-docker
    -just stop-syftbox-server
    @if [ -n "${CI:-}" ] && [ -d {{TEST_CLIENTS_DOCKER_DIR}} ]; then \
        echo "In CI environment, using sudo to remove root-owned files"; \
        sudo rm -rf {{TEST_CLIENTS_DOCKER_DIR}}; \
    else \
        rm -rf {{TEST_CLIENTS_DOCKER_DIR}}; \
    fi
    @echo "Cleanup complete!"

# ============= SBENV LOCAL CLIENT COMMANDS =============

# Build syftbox binary from source
build-syftbox-binary:
    #!/usr/bin/env bash
    set -e
    echo "Building SyftBox binary for current platform..."

    # Clean and create output directory
    rm -rf {{SYFTBOX_DIR}}/.out
    mkdir -p {{SYFTBOX_DIR}}/.out

    # Build directly with go build (blocking)
    echo "Building SyftBox binary ($(go env GOOS)/$(go env GOARCH))..."
    cd {{SYFTBOX_DIR}}

    # Get version info for build
    SYFTBOX_VERSION=$(git describe --tags --always 2>/dev/null || echo "dev")
    BUILD_COMMIT=$(git rev-parse --short HEAD)
    BUILD_DATE=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    BUILD_LD_FLAGS="-s -w -X github.com/openmined/syftbox/internal/version.Version=$SYFTBOX_VERSION -X github.com/openmined/syftbox/internal/version.Revision=$BUILD_COMMIT -X github.com/openmined/syftbox/internal/version.BuildDate=$BUILD_DATE"

    export GOOS=$(go env GOOS)
    export GOARCH=$(go env GOARCH)
    export CGO_ENABLED=0

    # Enable CGO for darwin
    if [ "$GOOS" = "darwin" ]; then
        export CGO_ENABLED=1
    fi

    # Run blocking go build (we're already in syftbox dir)
    go build -trimpath --tags="go_json nomsgpack" \
        -ldflags="$BUILD_LD_FLAGS" \
        -o .out/syftbox_client_${GOOS}_${GOARCH} ./cmd/client

    echo "Build completed. Looking for binary..."
    # List what was created (we're in syftbox dir, so use .out directly)
    ls -la .out/

    # Find the built binary (use .out since we're in syftbox dir)
    SYFTBOX_BIN=$(find .out -type f -name "syftbox_client_*" | head -1)

    if [ -n "$SYFTBOX_BIN" ]; then
        # Create absolute path for use by sbenv
        FULL_PATH="$(pwd)/$SYFTBOX_BIN"
        echo "Found SyftBox binary at: $FULL_PATH"

        # Test that it's executable and correct architecture
        if [ -x "$FULL_PATH" ]; then
            echo "Binary is executable, testing version..."
            "$FULL_PATH" --version || true
        fi

        # Export for use in subsequent commands (go back to workspace dir)
        cd ..
        echo "export SYFTBOX_BINARY_PATH=\"$FULL_PATH\"" > .syftbox_binary_path
        echo "SyftBox binary will be used by sbenv: $FULL_PATH"
    else
        echo "Error: Could not find SyftBox binary after build"
        exit 1
    fi

# Prepare sbenv (depends on syftbox binary being built)
prepare-sbenv: build-syftbox-binary
    @echo "sbenv will be run via cargo run shim for auto-rebuild"

# Verify sbenv environment for a client
verify-sbenv-client email:
    @echo "Verifying sbenv environment for {{email}}..."
    @if [ ! -d "{{TEST_CLIENTS_LOCAL_DIR}}/{{email}}" ]; then \
        echo "❌ Client directory does not exist: {{TEST_CLIENTS_LOCAL_DIR}}/{{email}}"; \
        exit 1; \
    fi
    @cd {{TEST_CLIENTS_LOCAL_DIR}}/{{email}} && \
    bash -c 'set +u; \
    echo "Activating sbenv environment..."; \
    eval "$(../../{{SBENV_DIR}}/sbenv activate --quiet)"; \
    echo ""; \
    echo "=== SYFTBOX Environment Variables ==="; \
    echo "All SYFTBOX_ variables:"; \
    env | grep "^SYFTBOX_" | sort || echo "No SYFTBOX_ variables found!"; \
    echo ""; \
    echo "=== Verification Results ==="; \
    ERRORS=0; \
    if [ -z "$SYFTBOX_BINARY" ]; then \
        echo "❌ SYFTBOX_BINARY is not set"; \
        ERRORS=$((ERRORS + 1)); \
    else \
        echo "✓ SYFTBOX_BINARY=$SYFTBOX_BINARY"; \
        if [ ! -f "$SYFTBOX_BINARY" ]; then \
            echo "  ❌ Binary does not exist at $SYFTBOX_BINARY"; \
            ERRORS=$((ERRORS + 1)); \
        else \
            echo "  ✓ Binary exists ($(du -h "$SYFTBOX_BINARY" | cut -f1))"; \
            if [ -x "$SYFTBOX_BINARY" ]; then \
                echo "  ✓ Binary is executable"; \
                $SYFTBOX_BINARY --version 2>&1 | head -1 || echo "  ❌ Could not run binary"; \
            else \
                echo "  ❌ Binary is not executable"; \
                ERRORS=$((ERRORS + 1)); \
            fi; \
        fi; \
    fi; \
    if [ -z "$SYFTBOX_VERSION" ]; then \
        echo "❌ SYFTBOX_VERSION is not set"; \
        ERRORS=$((ERRORS + 1)); \
    else \
        echo "✓ SYFTBOX_VERSION=$SYFTBOX_VERSION"; \
    fi; \
    if [ -z "$SYFTBOX_DATA_DIR" ]; then \
        echo "❌ SYFTBOX_DATA_DIR is not set"; \
        ERRORS=$((ERRORS + 1)); \
    else \
        echo "✓ SYFTBOX_DATA_DIR=$SYFTBOX_DATA_DIR"; \
    fi; \
    if [ -z "$SYFTBOX_EMAIL" ]; then \
        echo "❌ SYFTBOX_EMAIL is not set"; \
        ERRORS=$((ERRORS + 1)); \
    else \
        echo "✓ SYFTBOX_EMAIL=$SYFTBOX_EMAIL"; \
    fi; \
    if [ -z "$SYFTBOX_SERVER_URL" ]; then \
        echo "❌ SYFTBOX_SERVER_URL is not set"; \
        ERRORS=$((ERRORS + 1)); \
    else \
        echo "✓ SYFTBOX_SERVER_URL=$SYFTBOX_SERVER_URL"; \
    fi; \
    echo ""; \
    if [ $ERRORS -eq 0 ]; then \
        echo "✅ All environment variables are properly set!"; \
    else \
        echo "❌ Found $ERRORS errors in environment setup"; \
        exit 1; \
    fi'

# Setup local client environment with sbenv
setup-sbenv-client email: prepare-sbenv
    @echo "Setting up sbenv for {{email}}..."
    @mkdir -p {{TEST_CLIENTS_LOCAL_DIR}}/{{email}}
    @# Source the binary path if available
    @if [ -f .syftbox_binary_path ]; then \
        echo "Found .syftbox_binary_path file, contents:"; \
        cat .syftbox_binary_path; \
        source .syftbox_binary_path; \
        echo "After sourcing, SYFTBOX_BINARY_PATH=$SYFTBOX_BINARY_PATH"; \
        if [ -z "$SYFTBOX_BINARY_PATH" ]; then \
            echo "ERROR: SYFTBOX_BINARY_PATH is empty after sourcing!"; \
            exit 1; \
        fi; \
        echo "Using SyftBox binary: $SYFTBOX_BINARY_PATH"; \
        if [ ! -f "$SYFTBOX_BINARY_PATH" ]; then \
            echo "ERROR: Binary does not exist at $SYFTBOX_BINARY_PATH"; \
            exit 1; \
        fi; \
        echo "Binary exists, proceeding with sbenv init..."; \
        cd {{TEST_CLIENTS_LOCAL_DIR}}/{{email}} && \
        echo "Running command: ../../{{SBENV_DIR}}/sbenv init --dev --server-url http://localhost:8080 --email {{email}} --binary $SYFTBOX_BINARY_PATH"; \
        ../../{{SBENV_DIR}}/sbenv init --dev --server-url http://localhost:8080 --email {{email}} --binary "$SYFTBOX_BINARY_PATH"; \
        echo "sbenv init completed"; \
        echo "Checking .syftbox directory contents:"; \
        ls -la .syftbox/ 2>/dev/null || echo "No .syftbox directory yet"; \
        if [ -f .syftbox/config.json ]; then echo "Config:"; cat .syftbox/config.json; fi; \
    else \
        echo "No .syftbox_binary_path file found, running sbenv init without binary flag"; \
        cd {{TEST_CLIENTS_LOCAL_DIR}}/{{email}} && \
        ../../{{SBENV_DIR}}/sbenv init --dev --server-url http://localhost:8080 --email {{email}}; \
    fi

# Start local client with sbenv (uses environment from sbenv init)
start-sbenv-client email:
    @echo "Starting sbenv client {{email}}..."
    @cd {{TEST_CLIENTS_LOCAL_DIR}}/{{email}} && \
    bash -c 'set +u; \
    echo "Current directory: $(pwd)"; \
    echo "Activating sbenv environment..."; \
    eval "$(../../{{SBENV_DIR}}/sbenv activate --quiet)"; \
    echo ""; \
    echo "=== Verifying SYFTBOX Environment Variables ==="; \
    echo "All SYFTBOX_ variables:"; \
    env | grep "^SYFTBOX_" | sort || echo "No SYFTBOX_ variables found!"; \
    echo ""; \
    echo "Checking required variables:"; \
    MISSING_VARS=""; \
    if [ -z "$SYFTBOX_BINARY" ]; then \
        echo "❌ SYFTBOX_BINARY is not set"; \
        MISSING_VARS="yes"; \
    else \
        echo "✓ SYFTBOX_BINARY=$SYFTBOX_BINARY"; \
        if [ ! -f "$SYFTBOX_BINARY" ]; then \
            echo "  ⚠️  WARNING: Binary does not exist at $SYFTBOX_BINARY"; \
        else \
            echo "  ✓ Binary exists and is $(du -h "$SYFTBOX_BINARY" | cut -f1)"; \
        fi; \
    fi; \
    if [ -z "$SYFTBOX_VERSION" ]; then \
        echo "❌ SYFTBOX_VERSION is not set"; \
        MISSING_VARS="yes"; \
    else \
        echo "✓ SYFTBOX_VERSION=$SYFTBOX_VERSION"; \
    fi; \
    if [ -z "$SYFTBOX_DATA_DIR" ]; then \
        echo "❌ SYFTBOX_DATA_DIR is not set"; \
        MISSING_VARS="yes"; \
    else \
        echo "✓ SYFTBOX_DATA_DIR=$SYFTBOX_DATA_DIR"; \
        if [ ! -d "$SYFTBOX_DATA_DIR" ]; then \
            echo "  ⚠️  WARNING: Data dir does not exist at $SYFTBOX_DATA_DIR"; \
        fi; \
    fi; \
    if [ -z "$SYFTBOX_EMAIL" ]; then \
        echo "❌ SYFTBOX_EMAIL is not set"; \
        MISSING_VARS="yes"; \
    else \
        echo "✓ SYFTBOX_EMAIL=$SYFTBOX_EMAIL"; \
    fi; \
    if [ -z "$SYFTBOX_SERVER_URL" ]; then \
        echo "❌ SYFTBOX_SERVER_URL is not set"; \
        MISSING_VARS="yes"; \
    else \
        echo "✓ SYFTBOX_SERVER_URL=$SYFTBOX_SERVER_URL"; \
    fi; \
    echo ""; \
    if [ -n "$MISSING_VARS" ]; then \
        echo "⚠️  Critical environment variables are missing!"; \
        echo "Cannot proceed with starting client."; \
        exit 1; \
    fi; \
    echo "✓ All required environment variables are set"; \
    echo ""; \
    if [ -n "$SYFTBOX_BINARY" ] && [ -f "$SYFTBOX_BINARY" ]; then \
        echo "Creating temporary symlink for syftbox binary..."; \
        ln -sf "$SYFTBOX_BINARY" ./syftbox; \
        export PATH="$(pwd):$PATH"; \
        echo "Updated PATH to include: $(pwd)"; \
        echo "Syftbox in path: $(which syftbox 2>/dev/null || echo not found)"; \
    fi; \
    echo ""; \
    echo "Starting sbenv client..."; \
    ../../{{SBENV_DIR}}/sbenv start --skip-login-check 2>&1; \
    echo ""; \
    echo "Checking daemon log..."; \
    if [ -f .syftbox/daemon.log ]; then \
        echo "Last 20 lines of daemon.log:"; \
        tail -20 .syftbox/daemon.log; \
    else \
        echo "No daemon.log found"; \
    fi; \
    echo "sbenv start completed with code: $?"' &

# Setup and start sbenv client 1
start-syftbox-client1-local:
    @just setup-sbenv-client "{{TEST_CLIENT1_EMAIL}}"
    @just start-sbenv-client "{{TEST_CLIENT1_EMAIL}}"

# Setup and start sbenv client 2
start-syftbox-client2-local:
    @just setup-sbenv-client "{{TEST_CLIENT2_EMAIL}}"
    @just start-sbenv-client "{{TEST_CLIENT2_EMAIL}}"

# Setup and start sbenv bad client
start-syftbox-bad-client-local:
    @just setup-sbenv-client "{{TEST_BAD_CLIENT_EMAIL}}"
    @just start-sbenv-client "{{TEST_BAD_CLIENT_EMAIL}}"

# Start all local clients (including bad client)
start-syftbox-all-clients-local:
    @echo "Cleaning up old test data..."
    @rm -rf {{TEST_CLIENTS_LOCAL_DIR}}
    @echo "Starting client1..."
    @just start-syftbox-client1-local
    @sleep 5
    @echo "Starting client2..."
    @just start-syftbox-client2-local
    @sleep 5
    @echo "Starting bad client..."
    @just start-syftbox-bad-client-local
    @echo "Waiting for local clients to initialize..."
    @sleep 10
    @echo "Checking if clients are running..."
    @ps aux | grep -E "syftbox|sbenv" | grep -v grep || echo "No syftbox processes found"
    @echo "Checking client logs..."
    @for client in "{{TEST_CLIENT1_EMAIL}}" "{{TEST_CLIENT2_EMAIL}}" "{{TEST_BAD_CLIENT_EMAIL}}"; do \
        echo "=== Logs for $$client ==="; \
        if [ -f "{{TEST_CLIENTS_LOCAL_DIR}}/$$client/.syftbox/daemon.log" ]; then \
            tail -10 "{{TEST_CLIENTS_LOCAL_DIR}}/$$client/.syftbox/daemon.log"; \
        else \
            echo "No daemon.log found"; \
        fi; \
    done

# Start both local clients (original, without bad client)
start-syftbox-clients-local:
    @echo "Cleaning up old test data..."
    @rm -rf {{TEST_CLIENTS_LOCAL_DIR}}
    @just start-syftbox-client1-local
    @sleep 3
    @just start-syftbox-client2-local
    @echo "Waiting for local clients to initialize..."
    @sleep 10
    @echo "Checking if clients are running..."
    @ps aux | grep -E "syftbox|sbenv" | grep -v grep || echo "No syftbox processes found"

# Stop all local clients
stop-syftbox-all-clients-local:
    @echo "Stopping all local sbenv clients..."
    -cd {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_CLIENT1_EMAIL}} && \
        ../../{{SBENV_DIR}}/sbenv stop 2>/dev/null || true
    -cd {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_CLIENT2_EMAIL}} && \
        ../../{{SBENV_DIR}}/sbenv stop 2>/dev/null || true
    -cd {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_BAD_CLIENT_EMAIL}} && \
        ../../{{SBENV_DIR}}/sbenv stop 2>/dev/null || true

# Stop local clients (original)
stop-syftbox-clients-local:
    @echo "Stopping local sbenv clients..."
    -cd {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_CLIENT1_EMAIL}} && \
        ../../{{SBENV_DIR}}/sbenv stop 2>/dev/null || true
    -cd {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_CLIENT2_EMAIL}} && \
        ../../{{SBENV_DIR}}/sbenv stop 2>/dev/null || true

# Run all integration tests with local clients
run-integration-tests-local:
    @echo "Building BioVault binary..."
    cd cli && cargo build --release
    @echo "Running integration tests with local clients..."
    cd cli && \
    SYFTBOX_SERVER_URL=http://localhost:8080 \
    SYFTBOX_CLIENT1_EMAIL={{TEST_CLIENT1_EMAIL}} \
    SYFTBOX_CLIENT2_EMAIL={{TEST_CLIENT2_EMAIL}} \
    SYFTBOX_BAD_CLIENT_EMAIL={{TEST_BAD_CLIENT_EMAIL}} \
    TEST_CLIENTS_DIR=../{{TEST_CLIENTS_LOCAL_DIR}} \
    TEST_MODE=local \
    cargo test --features e2e-tests --test '*' -- --ignored --nocapture

# Full integration test with local clients
test-integration-local cleanup="true" reset="false":
    @echo "=== Starting Integration Test with Local Clients ==="
    @if [ "{{reset}}" = "true" ] || [ "{{reset}}" = "--reset" ]; then \
        echo "Reset requested, cleaning data..."; \
        just quick-reset; \
    fi
    @echo "1. Checking /etc/hosts for minio..."
    @if ! grep -q "^[^#]*127.0.0.1.*minio" /etc/hosts; then \
        echo "❌ minio not found in /etc/hosts. Running: just add-minio-host"; \
        just add-minio-host; \
    else \
        echo "✓ minio found in /etc/hosts"; \
    fi
    @echo ""
    @echo "2. Starting SyftBox server..."
    @just start-syftbox-server
    @echo ""
    @echo "3. Waiting for server..."
    @just wait-for-server
    @echo ""
    @echo "4. Starting local SyftBox clients with sbenv (including bad client)..."
    @just start-syftbox-all-clients-local
    @echo ""
    @echo "5. Running integration tests..."
    @if just run-integration-tests-local; then \
        echo ""; \
        if [ "{{cleanup}}" = "true" ]; then \
            echo "✅ All tests passed! Cleaning up..."; \
            just stop-syftbox-all-clients-local; \
            just stop-syftbox-server; \
            rm -rf {{TEST_CLIENTS_LOCAL_DIR}}; \
            echo "✅ Cleanup complete!"; \
        else \
            echo "✅ All tests passed! Leaving services running for inspection."; \
            echo ""; \
            echo "Services are still running. You can:"; \
            echo "  - Check server logs: docker logs syftbox-server"; \
            echo "  - Check client1 logs: cat {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_CLIENT1_EMAIL}}/.syftbox/logs/client.log"; \
            echo "  - Check client2 logs: cat {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_CLIENT2_EMAIL}}/.syftbox/logs/client.log"; \
            echo "  - Check bad client logs: cat {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_BAD_CLIENT_EMAIL}}/.syftbox/logs/client.log"; \
            echo "  - Explore files: ls -la {{TEST_CLIENTS_LOCAL_DIR}}/"; \
            echo ""; \
            echo "When done, run:"; \
            echo "  just cleanup-integration-test-local"; \
        fi \
    else \
        echo ""; \
        echo "❌ Tests failed! Leaving services running for debugging."; \
        echo ""; \
        echo "Debug commands:"; \
        echo "  - Server logs: docker logs syftbox-server"; \
        echo "  - Client1 logs: cat {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_CLIENT1_EMAIL}}/.syftbox/logs/client.log"; \
        echo "  - Client2 logs: cat {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_CLIENT2_EMAIL}}/.syftbox/logs/client.log"; \
        echo "  - Bad client logs: cat {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_BAD_CLIENT_EMAIL}}/.syftbox/logs/client.log"; \
        echo "  - Check files: ls -la {{TEST_CLIENTS_LOCAL_DIR}}/"; \
        echo ""; \
        echo "When done debugging, run:"; \
        echo "  just cleanup-integration-test-local"; \
        exit 1; \
    fi

# Run local integration test without cleanup
test-integration-local-inspect:
    just test-integration-local false false

# Run local integration test with reset
test-integration-local-reset:
    just test-integration-local true --reset

# Clean up local integration test resources
cleanup-integration-test-local:
    @echo "Cleaning up local integration test resources..."
    -just stop-syftbox-all-clients-local
    -just stop-syftbox-server
    -rm -rf {{TEST_CLIENTS_LOCAL_DIR}}
    @echo "Cleanup complete!"


# Show status of services
status:
    @echo "=== Service Status ==="
    @printf "Server: "
    @if [ "$(just server-is-running)" = "true" ]; then echo "✓ Running"; else echo "✗ Stopped"; fi
    @printf "MinIO: "
    @if [ "$(just minio-is-running)" = "true" ]; then echo "✓ Running"; else echo "✗ Stopped"; fi
    @printf "Docker Clients: "
    @docker ps --format 'table {{{{.Names}}}}' | grep -c syftbox-client || echo "0"
    @echo ""
    @echo "Test Directories:"
    @if [ -d {{TEST_CLIENTS_DOCKER_DIR}} ]; then echo "  ✓ {{TEST_CLIENTS_DOCKER_DIR}} exists"; else echo "  ✗ {{TEST_CLIENTS_DOCKER_DIR}} missing"; fi
    @if [ -d {{TEST_CLIENTS_LOCAL_DIR}} ]; then echo "  ✓ {{TEST_CLIENTS_LOCAL_DIR}} exists"; else echo "  ✗ {{TEST_CLIENTS_LOCAL_DIR}} missing"; fi
