# COLOC-flow

A Nextflow pipeline for colocalization analysis to identify shared genetic variants between GWAS and eQTL datasets.

## Overview

This pipeline implements colocalization analysis using the `coloc` R package to identify genetic variants that influence both disease risk and gene expression. The workflow processes GWAS summary statistics and eQTL data, performs colocalization analysis, and generates visualizations and reports.

## Requirements

- Nextflow (>=20.10.0)
- Singularity or Docker
- R (with coloc, data.table, and ggplot2 packages)

## Usage

```bash
nextflow run main.nf --gwas_file <path_to_gwas_file> --eqtl_file <path_to_eqtl_file>
```

## Profiles

The pipeline comes with predefined configuration profiles:

- `standard`: Default profile that runs locally
- `offline`: Uses Docker containers and runs locally
- `cluster`: Configured for PBS Pro job scheduler with Singularity containers

Run with a specific profile:

```bash
nextflow run main.nf -profile cluster --gwas_file <path> --eqtl_file <path>
```

## Pipeline Steps

1. **PREPARE_GWAS**: Process and format GWAS summary statistics
2. **PROCESS_EQTL**: Process and format eQTL data
3. **RUN_COLOC**: Perform colocalization analysis using the coloc R package

