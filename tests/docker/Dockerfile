# Emulate Google Colab environment for testing
FROM ubuntu:22.04

# Set environment variables to emulate Colab
ENV COLAB_RELEASE_TAG=release-colab_20250822-060052_RC00
ENV DEBIAN_FRONTEND=noninteractive

# Update and install basic packages
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    ca-certificates \
    gnupg \
    lsb-release \
    build-essential \
    pkg-config \
    libssl-dev \
    git \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Install Java 11 (default in Colab)
RUN apt-get update && apt-get install -y \
    openjdk-11-jdk-headless \
    && rm -rf /var/lib/apt/lists/*

# Install Rust toolchain for building the bv binary
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install uv (Python package manager) to match real Colab environment
# Real Colab has uv 0.9.0 pre-installed at /usr/local/bin/uv
RUN curl -LsSf https://astral.sh/uv/install.sh | sh && \
    mv /root/.cargo/bin/uv /usr/local/bin/uv && \
    chmod +x /usr/local/bin/uv && \
    /usr/local/bin/uv --version

# Ensure /usr/local/bin is in PATH for all processes (matches real Colab)
ENV PATH="/usr/local/bin:${PATH}"

# Create working directory
WORKDIR /workspace

# Create a directory for the built binary
RUN mkdir -p /usr/local/bin

# Set up entrypoint
ENTRYPOINT ["/bin/bash"]