FROM rocker/r-ver:4.3.2

LABEL maintainer="Reuben Duncan <reuben.duncan25@outlook.com>"
LABEL description="Alpha diversity analysis — R environment"

# System dependencies for R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libgit2-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libhdf5-dev \
        libzstd-dev \
        liblz4-dev \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN Rscript -e "\
    install.packages( \
        c('optparse', 'stringr', 'data.table', 'vegan', 'ape', 'car'), \
        repos = 'https://cloud.r-project.org', \
        Ncpus = max(1L, parallel::detectCores() - 1L) \
    )"

# Install Bioconductor packages (biomformat for BIOM I/O; picante for Faith's PD)
RUN Rscript -e "\
    if (!requireNamespace('BiocManager', quietly = TRUE)) \
        install.packages('BiocManager', repos = 'https://cloud.r-project.org'); \
    BiocManager::install(c('biomformat', 'picante'), ask = FALSE, update = FALSE)"

# Install arrow (pre-built C++ library; LIBARROW_BINARY avoids 30-min source compile)
RUN LIBARROW_BINARY=true Rscript -e "\
    install.packages('arrow', repos='https://cloud.r-project.org', \
        Ncpus=max(1L, parallel::detectCores()-1L))"

# Copy R helper scripts into the container
COPY src/ /opt/ecology-scripts/src/

# Verify all required packages load cleanly
RUN Rscript -e "\
    library(optparse); library(stringr); library(data.table); \
    library(vegan); library(biomformat); library(ape); library(picante); \
    library(car); library(arrow); \
    message('All packages loaded successfully.')"

WORKDIR /data
