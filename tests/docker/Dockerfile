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

# Create working directory
WORKDIR /workspace

# Create a directory for the built binary
RUN mkdir -p /usr/local/bin

# Set up entrypoint
ENTRYPOINT ["/bin/bash"]