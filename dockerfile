FROM rocker/tidyverse:latest

LABEL maintainer="Johnson Lab <your.email@example.com>"
LABEL description="Docker image for CAUSAL-flow pipeline"

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
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install plink directly
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip \
    && unzip plink_linux_x86_64_20241022.zip -d /usr/local/bin/ \
    && rm plink_linux_x86_64_20241022.zip \
    && chmod +x /usr/local/bin/plink

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
RUN R -e "devtools::install_github('MRCIEU/genetics.binaRies')"

RUN R -e "install.packages('MendelianRandomization')"
# Create working directory
RUN mkdir -p /app/CAUSAL-flow

#get reference data 
RUN wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz \
    && tar -xzf 1kg.v3.tgz -C /app/CAUSAL-flow/ \
    && rm 1kg.v3.tgz

WORKDIR /app/CAUSAL-flow

# Set the default command
# This allows overriding with "docker run -it causal-flow:latest bash"
CMD ["Rscript", "--help"]