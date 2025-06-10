# CAUSAL-flow

A Nextflow pipeline for colocalization analysis to identify shared genetic variants between GWAS and eQTL datasets.

## Overview

This pipeline implements colocalization analysis using the `coloc` R package to identify genetic variants that influence both disease risk and gene expression. The workflow processes GWAS summary statistics and eQTL data, performs colocalization analysis, and generates visualizations and reports.

## Requirements

- Nextflow (>=20.10.0)
- Singularity or Docker
- R (with coloc, data.table, and ggplot2 packages)

## Usage

```bash
nextflow run main.nf \
  --gwas_file <path/to/gwas_summary.txt> \
  --gwas_name <name_of_gwas_study> \
  --eqtl_data <path/to/eqtl_data.rds> \
  --gene_location_file <path/to/gene_locations.csv> \
  --outdir <output_directory>
```

### Input Files

1. **GWAS Summary Statistics**:
   - Must contain columns: SNP, CHR, BP, P, BETA, SE

2. **eQTL Data**:
   - RDS file containing eQTL data with columns: SNP, gene, beta, t.stat

3. **Gene Location File**:
   - CSV file containing gene coordinates with columns: gene_id, chr, start, end

### Example Command

```bash
nextflow run main.nf \
  --gwas_file /path/to/Epilepsy_GWAS.txt \
  --gwas_name Epilepsy \
  --eqtl_data /path/to/Brain_eQTL_data.rds \
  --gene_location_file /path/to/gene_annotations_hg38.csv \
  --outdir results/epilepsy_coloc
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

