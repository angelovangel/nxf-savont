# Build stage
FROM rust:bookworm AS builder
ARG TARGETARCH

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    cmake \
    libclang-dev \
    pkg-config \
    wget \
    tar \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/bluenote-1577/savont && \
    cd savont && \
    cargo install --path . && \
    cd .. && \
    git clone https://github.com/angelovangel/faster && \
    cd faster && \
    cargo build --release

# Taxonkit installation
RUN ARCH=${TARGETARCH:-$(uname -m)} && \
    case "$ARCH" in \
    arm64|aarch64) TARGET="arm64" ;; \
    amd64|x86_64) TARGET="amd64" ;; \
    *) echo "Unsupported architecture: $ARCH" && exit 1 ;; \
    esac && \
    wget https://github.com/shenwei356/taxonkit/releases/download/v0.20.0/taxonkit_linux_${TARGET}.tar.gz -O taxonkit.tar.gz && \
    tar -zxvf taxonkit.tar.gz taxonkit && \
    mv taxonkit /usr/local/cargo/bin/ && \
    rm taxonkit.tar.gz

# NCBI Taxonomy database (for taxonkit)
RUN wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && \
    mkdir -p /opt/taxonkit && \
    tar -zxvf taxdump.tar.gz -C /opt/taxonkit names.dmp nodes.dmp delnodes.dmp merged.dmp && \
    rm taxdump.tar.gz

# Savont databases
RUN mkdir -p /databases && \
    /usr/local/cargo/bin/savont download --location /databases --emu-db && \
    /usr/local/cargo/bin/savont download --location /databases --silva-db

# Final stage
FROM debian:bookworm-slim
LABEL name="aangeloo/nxf-savont"
LABEL maintainer="https://github.com/angelovangel"

COPY --from=builder /usr/local/cargo/bin/savont /usr/local/bin/savont
COPY --from=builder /faster/target/release/faster /usr/local/bin/faster
COPY --from=builder /usr/local/cargo/bin/taxonkit /usr/local/bin/taxonkit
COPY --from=builder /opt/taxonkit /opt/taxonkit
COPY --from=builder /databases /databases

# 
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/* \
    && chmod -R a+rx /opt/taxonkit
