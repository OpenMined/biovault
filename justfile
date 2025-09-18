DOCKER_DIR := "./syftbox/docker"
SYFTBOX_DIR := "./syftbox"
SBENV_DIR := "./sbenv"
TEST_CLIENTS_DOCKER_DIR := "./test-clients-docker"
TEST_CLIENTS_LOCAL_DIR := "./test-clients-local"

# Test client configurations
TEST_CLIENT1_EMAIL := "client1@syftbox.net"
TEST_CLIENT2_EMAIL := "client2@syftbox.net"
TEST_CLIENT1_NAME := "client1-syftbox-net"
TEST_CLIENT2_NAME := "client2-syftbox-net"
TEST_CLIENT1_PORT := "7938"
TEST_CLIENT2_PORT := "7939"

start-syftbox-server:
    @echo "Starting SyftBox server with MinIO..."
    cd {{DOCKER_DIR}} && COMPOSE_BAKE=true docker compose up -d --build minio server
    @echo "Waiting for server to be ready..."
    @sleep 5
    @echo "Server started at http://localhost:8080"

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
    @chmod -R 777 {{TEST_CLIENTS_DOCKER_DIR}}
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

# Start both docker clients
start-syftbox-clients-docker:
    @echo "Cleaning up old test data..."
    @rm -rf {{TEST_CLIENTS_DOCKER_DIR}}
    @just start-syftbox-client1-docker
    @sleep 2
    @just start-syftbox-client2-docker

# Stop docker clients
stop-syftbox-clients-docker:
    -docker stop syftbox-client-{{TEST_CLIENT1_NAME}}
    -docker stop syftbox-client-{{TEST_CLIENT2_NAME}}
    -docker rm syftbox-client-{{TEST_CLIENT1_NAME}}
    -docker rm syftbox-client-{{TEST_CLIENT2_NAME}}

# Wait for server to be ready
wait-for-server:
    @echo "Waiting for server to be ready..."
    @timeout 60 bash -c 'until curl -f http://localhost:8080 2>/dev/null; do echo "Waiting..."; sleep 2; done'
    @echo "Server is ready!"

# Run integration tests with docker clients
run-integration-tests-docker:
    @echo "Running integration tests with docker clients..."
    cd cli && \
    SYFTBOX_SERVER_URL=http://localhost:8080 \
    SYFTBOX_CLIENT1_EMAIL={{TEST_CLIENT1_EMAIL}} \
    SYFTBOX_CLIENT2_EMAIL={{TEST_CLIENT2_EMAIL}} \
    TEST_CLIENTS_DIR=../{{TEST_CLIENTS_DOCKER_DIR}} \
    TEST_MODE=docker \
    cargo test --test '*' -- --ignored --nocapture

# Full integration test with Docker - starts everything, runs tests, cleans up on success
test-integration-docker cleanup="true":
    @echo "=== Starting Integration Test with Docker ==="
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
    @echo "4. Starting SyftBox docker clients..."
    @just start-syftbox-clients-docker
    @echo ""
    @echo "5. Waiting for clients to initialize..."
    @sleep 5
    @echo ""
    @echo "6. Running integration tests..."
    @if just run-integration-tests-docker; then \
        echo ""; \
        if [ "{{cleanup}}" = "true" ]; then \
            echo "✅ All tests passed! Cleaning up..."; \
            just stop-syftbox-clients-docker; \
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
        echo "  - Check files: ls -la {{TEST_CLIENTS_DOCKER_DIR}}/"; \
        echo ""; \
        echo "When done debugging, run:"; \
        echo "  just cleanup-integration-test"; \
        exit 1; \
    fi

# Run integration test without cleanup (keeps services running after)
test-integration-docker-inspect:
    just test-integration-docker false

# Clean up docker integration test resources
cleanup-integration-test-docker:
    @echo "Cleaning up docker integration test resources..."
    -just stop-syftbox-clients-docker
    -just stop-syftbox-server
    -rm -rf {{TEST_CLIENTS_DOCKER_DIR}}
    @echo "Cleanup complete!"

# ============= SBENV LOCAL CLIENT COMMANDS =============

# Build sbenv tool
build-sbenv:
    @echo "Building sbenv..."
    @if [ ! -f "{{SBENV_DIR}}/cli/target/release/sbenv" ]; then \
        cd {{SBENV_DIR}}/cli && cargo build --release; \
    else \
        echo "sbenv already built"; \
    fi

# Setup local client environment with sbenv
setup-sbenv-client email:
    @echo "Setting up sbenv for {{email}}..."
    @mkdir -p {{TEST_CLIENTS_LOCAL_DIR}}/{{email}}
    cd {{TEST_CLIENTS_LOCAL_DIR}}/{{email}} && \
    ../../{{SBENV_DIR}}/cli/target/release/sbenv init --dev --server-url http://localhost:8080 --email {{email}}

# Start local client with sbenv (uses environment from sbenv init)
start-sbenv-client email:
    @echo "Starting sbenv client {{email}}..."
    cd {{TEST_CLIENTS_LOCAL_DIR}}/{{email}} && \
    set +u && \
    eval "$(../../{{SBENV_DIR}}/cli/target/release/sbenv activate --quiet)" && \
    ../../{{SBENV_DIR}}/cli/target/release/sbenv start &

# Setup and start sbenv client 1
start-syftbox-client1-local: build-sbenv
    @just setup-sbenv-client "{{TEST_CLIENT1_EMAIL}}"
    @just start-sbenv-client "{{TEST_CLIENT1_EMAIL}}"

# Setup and start sbenv client 2
start-syftbox-client2-local: build-sbenv
    @just setup-sbenv-client "{{TEST_CLIENT2_EMAIL}}"
    @just start-sbenv-client "{{TEST_CLIENT2_EMAIL}}"

# Start both local clients
start-syftbox-clients-local: build-sbenv
    @echo "Cleaning up old test data..."
    @rm -rf {{TEST_CLIENTS_LOCAL_DIR}}
    @just start-syftbox-client1-local
    @sleep 2
    @just start-syftbox-client2-local
    @echo "Waiting for local clients to initialize..."
    @sleep 5

# Stop local clients
stop-syftbox-clients-local:
    @echo "Stopping local sbenv clients..."
    -cd {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_CLIENT1_EMAIL}} && \
        ../../{{SBENV_DIR}}/cli/target/release/sbenv stop 2>/dev/null || true
    -cd {{TEST_CLIENTS_LOCAL_DIR}}/{{TEST_CLIENT2_EMAIL}} && \
        ../../{{SBENV_DIR}}/cli/target/release/sbenv stop 2>/dev/null || true

# Run integration tests with local clients
run-integration-tests-local:
    @echo "Running integration tests with local clients..."
    cd cli && \
    SYFTBOX_SERVER_URL=http://localhost:8080 \
    SYFTBOX_CLIENT1_EMAIL={{TEST_CLIENT1_EMAIL}} \
    SYFTBOX_CLIENT2_EMAIL={{TEST_CLIENT2_EMAIL}} \
    TEST_CLIENTS_DIR=../{{TEST_CLIENTS_LOCAL_DIR}} \
    TEST_MODE=local \
    cargo test --test '*' -- --ignored --nocapture

# Full integration test with local clients
test-integration-local cleanup="true":
    @echo "=== Starting Integration Test with Local Clients ==="
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
    @echo "4. Starting local SyftBox clients with sbenv..."
    @just start-syftbox-clients-local
    @echo ""
    @echo "5. Running integration tests..."
    @if just run-integration-tests-local; then \
        echo ""; \
        if [ "{{cleanup}}" = "true" ]; then \
            echo "✅ All tests passed! Cleaning up..."; \
            just stop-syftbox-clients-local; \
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
        echo "  - Check files: ls -la {{TEST_CLIENTS_LOCAL_DIR}}/"; \
        echo ""; \
        echo "When done debugging, run:"; \
        echo "  just cleanup-integration-test-local"; \
        exit 1; \
    fi

# Run local integration test without cleanup
test-integration-local-inspect:
    just test-integration-local false

# Clean up local integration test resources
cleanup-integration-test-local:
    @echo "Cleaning up local integration test resources..."
    -just stop-syftbox-clients-local
    -just stop-syftbox-server
    -rm -rf {{TEST_CLIENTS_LOCAL_DIR}}
    @echo "Cleanup complete!"
