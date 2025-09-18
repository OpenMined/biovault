DOCKER_DIR := "./syftbox/docker"
SYFTBOX_DIR := "./syftbox"
TEST_CLIENTS_DIR := "./test-clients"

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

# Generic client start command
start-syftbox-client email name port:
    @echo "Starting client {{email}}..."
    @mkdir -p {{TEST_CLIENTS_DIR}}
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
        -v "$(pwd)/{{TEST_CLIENTS_DIR}}:/data/clients" \
        syftbox-client {{email}}
    @echo "Client started at http://localhost:{{port}}"

# Start syftbox client 1
start-syftbox-client1:
    just start-syftbox-client "{{TEST_CLIENT1_EMAIL}}" "{{TEST_CLIENT1_NAME}}" "{{TEST_CLIENT1_PORT}}"

# Start syftbox client 2
start-syftbox-client2:
    just start-syftbox-client "{{TEST_CLIENT2_EMAIL}}" "{{TEST_CLIENT2_NAME}}" "{{TEST_CLIENT2_PORT}}"

# Start both syftbox clients
start-syftbox-clients:
    @echo "Cleaning up old test data..."
    @rm -rf {{TEST_CLIENTS_DIR}}
    @just start-syftbox-client1
    @sleep 2
    @just start-syftbox-client2

# Stop syftbox clients
stop-syftbox-clients:
    -docker stop syftbox-client-{{TEST_CLIENT1_NAME}}
    -docker stop syftbox-client-{{TEST_CLIENT2_NAME}}
    -docker rm syftbox-client-{{TEST_CLIENT1_NAME}}
    -docker rm syftbox-client-{{TEST_CLIENT2_NAME}}

# Wait for server to be ready
wait-for-server:
    @echo "Waiting for server to be ready..."
    @timeout 60 bash -c 'until curl -f http://localhost:8080 2>/dev/null; do echo "Waiting..."; sleep 2; done'
    @echo "Server is ready!"

# Run integration tests only (assumes services are running)
run-integration-tests:
    @echo "Running integration tests..."
    cd cli && \
    SYFTBOX_SERVER_URL=http://localhost:8080 \
    SYFTBOX_CLIENT1_EMAIL={{TEST_CLIENT1_EMAIL}} \
    SYFTBOX_CLIENT2_EMAIL={{TEST_CLIENT2_EMAIL}} \
    TEST_CLIENTS_DIR=../{{TEST_CLIENTS_DIR}} \
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
    @echo "4. Starting SyftBox clients..."
    @just start-syftbox-clients
    @echo ""
    @echo "5. Waiting for clients to initialize..."
    @sleep 5
    @echo ""
    @echo "6. Running integration tests..."
    @if just run-integration-tests; then \
        echo ""; \
        if [ "{{cleanup}}" = "true" ]; then \
            echo "✅ All tests passed! Cleaning up..."; \
            just stop-syftbox-clients; \
            just stop-syftbox-server; \
            rm -rf {{TEST_CLIENTS_DIR}}; \
            echo "✅ Cleanup complete!"; \
        else \
            echo "✅ All tests passed! Leaving services running for inspection."; \
            echo ""; \
            echo "Services are still running. You can:"; \
            echo "  - Check server logs: docker logs syftbox-server"; \
            echo "  - Check client1 logs: docker logs syftbox-client-{{TEST_CLIENT1_NAME}}"; \
            echo "  - Check client2 logs: docker logs syftbox-client-{{TEST_CLIENT2_NAME}}"; \
            echo "  - Explore files: ls -la {{TEST_CLIENTS_DIR}}/"; \
            echo ""; \
            echo "When done, run:"; \
            echo "  just cleanup-integration-test"; \
        fi \
    else \
        echo ""; \
        echo "❌ Tests failed! Leaving services running for debugging."; \
        echo ""; \
        echo "Debug commands:"; \
        echo "  - Server logs: docker logs syftbox-server"; \
        echo "  - Client1 logs: docker logs syftbox-client-{{TEST_CLIENT1_NAME}}"; \
        echo "  - Client2 logs: docker logs syftbox-client-{{TEST_CLIENT2_NAME}}"; \
        echo "  - Check files: ls -la {{TEST_CLIENTS_DIR}}/"; \
        echo ""; \
        echo "When done debugging, run:"; \
        echo "  just cleanup-integration-test"; \
        exit 1; \
    fi

# Run integration test without cleanup (keeps services running after)
test-integration-docker-inspect:
    just test-integration-docker false

# Clean up integration test resources
cleanup-integration-test:
    @echo "Cleaning up integration test resources..."
    -just stop-syftbox-clients
    -just stop-syftbox-server
    -rm -rf {{TEST_CLIENTS_DIR}}
    @echo "Cleanup complete!"
