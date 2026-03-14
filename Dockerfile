# Build stage
FROM rust:bookworm AS builder
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    cmake \
    libclang-dev \
    pkg-config \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/bluenote-1577/savont && \
    cd savont && \
    cargo install --path . && \
    cd .. && \
    git clone https://github.com/angelovangel/faster && \
    cd faster && \
    cargo build --release

RUN mkdir -p /databases && \
    /usr/local/cargo/bin/savont download --location /databases --emu-db && \
    /usr/local/cargo/bin/savont download --location /databases --silva-db

# Final stage
FROM debian:bookworm-slim
LABEL name="aangeloo/nxf-savont"
LABEL maintainer="https://github.com/angelovangel"


COPY --from=builder /usr/local/cargo/bin/savont /usr/local/bin/savont
COPY --from=builder /faster/target/release/faster /usr/local/bin/faster
COPY --from=builder /databases /databases

# 
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# 



