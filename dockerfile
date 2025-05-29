FROM rocker/tidyverse:latest

LABEL maintainer="Your Name <your.email@example.com>"
LABEL description="Docker image for COLOC-flow pipeline"

# Install system dependencies required for R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    zlib1g-dev \
    libgsl-dev \
    libsodium-dev \
    libgit2-dev \
    plink \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cran.r-project.org/')"

# Install required R packages from CRAN and Bioconductor
RUN R -e "install.packages(c( \
    'argparse', \
    'coloc', \
    'data.table', \
    'devtools', \
    'dplyr', \
    'ggplot2', \
    'ggrepel', \
    'gridExtra', \
    'hyprcoloc', \
    'markdown', \
    'MendelianRandomization', \
    'patchwork', \
    'plotly', \
    'rmarkdown', \
    'signif', \
    'optparse', \
    'readr', \
    'writexl' \
    ), repos='https://cran.r-project.org/')"

# Install ieugwasr
RUN R -e "devtools::install_github('MRCIEU/ieugwasr')"

# Install additional Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'GenomicRanges', \
    'rtracklayer', \
    'MatrixEQTL', \
    'snpStats', \
    'biomaRt', \
    'limma', \
    'edgeR' \
    ), ask=FALSE)"

# Install genetics.binaRies which is used in the cellCOLOC source
RUN R -e "devtools::install_github('explodecomputer/genetics.binaRies')"

# Create working directory
RUN mkdir -p /app/COLOC-flow
WORKDIR /app/COLOC-flow

# Add a default entrypoint
ENTRYPOINT ["Rscript"]

# Default command
CMD ["--help"]