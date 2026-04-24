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
        c('optparse', 'stringr', 'data.table', 'vegan'), \
        repos = 'https://cloud.r-project.org', \
        Ncpus = max(1L, parallel::detectCores() - 1L) \
    )"

# Install Bioconductor + phyloseq
RUN Rscript -e "\
    if (!requireNamespace('BiocManager', quietly = TRUE)) \
        install.packages('BiocManager', repos = 'https://cloud.r-project.org'); \
    BiocManager::install('phyloseq', ask = FALSE, update = FALSE)"

# Install arrow (pre-built C++ library; LIBARROW_BINARY avoids 30-min source compile)
RUN LIBARROW_BINARY=true Rscript -e "\
    install.packages('arrow', repos='https://cloud.r-project.org', \
        Ncpus=max(1L, parallel::detectCores()-1L))"

# Copy R helper scripts into the container
COPY src/ /opt/ecology-scripts/src/

# Verify all required packages load cleanly
RUN Rscript -e "\
    library(optparse); library(stringr); library(data.table); \
    library(vegan); library(phyloseq); library(arrow); \
    message('All packages loaded successfully.')"

WORKDIR /data
