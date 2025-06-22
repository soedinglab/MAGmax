FROM rust:1.78 as builder

WORKDIR /app

COPY Cargo.toml Cargo.lock ./
RUN mkdir src && echo 'fn main() {}' > src/main.rs && cargo build --release
RUN rm -rf src

COPY src ./src
RUN cargo build --release

FROM debian:bookworm-slim

# Micromamba setup
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH

# Install micromamba and dependencies
RUN apt-get update && apt-get install -y curl bzip2 ca-certificates && \
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | \
    tar -xvj -C /usr/local/bin/ --strip-components=1 bin/micromamba && \
    mkdir -p $MAMBA_ROOT_PREFIX

RUN micromamba install -y -n base -c bioconda -c conda-forge \
    skani spades seqtk && \
    micromamba clean -a -y

WORKDIR /app
COPY --from=builder /app/target/release/magmax /usr/local/bin/magmax

ENTRYPOINT ["magmax"]