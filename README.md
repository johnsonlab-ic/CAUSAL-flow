# CAUSAL-flow

A Nextflow pipeline for colocalization analysis to identify shared genetic variants between GWAS and eQTL datasets.

## Overview

This pipeline implements colocalization analysis using the `coloc` R package to identify genetic variants that influence both disease risk and gene expression. The workflow processes GWAS summary statistics and eQTL data, performs colocalization analysis, and generates visualizations and reports.

## Requirements

- Nextflow (>=20.10.0)
- Singularity or Docker
- R (with coloc, data.table, and ggplot2 packages)

## Usage

You can run the pipeline directly from GitHub:

```bash
nextflow run johnsonlab-ic/CAUSAL-flow \
  --gwas_file <path/to/gwas_summary.txt> \
  --gwas_name <name_of_gwas_study> \
  --eqtl_file <path/to/eqtl_data.rds> \
  --gene_location_file <path/to/gene_locations.csv> \
  --outdir <output_directory>
```

Or if you've cloned the repository:

```bash
nextflow run main.nf \
  --gwas_file <path/to/gwas_summary.txt> \
  --gwas_name <name_of_gwas_study> \
  --eqtl_file <path/to/eqtl_data.rds> \
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

### Example Commands

#### Single GWAS + Single eQTL
```bash
nextflow run johnsonlab-ic/CAUSAL-flow \
  --gwas_file /path/to/Epilepsy_GWAS.txt \
  --gwas_name Epilepsy \
  --eqtl_file /path/to/Neuron_eQTL_data.rds \
  --gene_location_file /path/to/gene_annotations_hg38.csv \
  --outdir results/epilepsy_coloc
```

#### Multiple GWAS + Single eQTL
```bash
nextflow run johnsonlab-ic/CAUSAL-flow \
  --gwas_dir /path/to/gwas_directory/ \
  --gwas_pattern "*.txt" \
  --eqtl_file /path/to/Neuron_eQTL_data.rds \
  --gene_location_file /path/to/gene_annotations_hg38.csv \
  --outdir results/multiple_gwas_coloc
```

#### Single GWAS + Multiple eQTLs
```bash
nextflow run johnsonlab-ic/CAUSAL-flow \
  --gwas_file /path/to/Epilepsy_GWAS.txt \
  --gwas_name Epilepsy \
  --eqtl_dir /path/to/eqtl_directory/ \
  --eqtl_pattern "*_cis_MatrixEQTLout.rds" \
  --gene_location_file /path/to/gene_annotations_hg38.csv \
  --outdir results/multiple_eqtl_coloc
```

#### Multiple GWAS + Multiple eQTLs
```bash
nextflow run johnsonlab-ic/CAUSAL-flow \
  --gwas_dir /path/to/gwas_directory/ \
  --gwas_pattern "*.txt" \
  --eqtl_dir /path/to/eqtl_directory/ \
  --eqtl_pattern "*_cis_MatrixEQTLout.rds" \
  --gene_location_file /path/to/gene_annotations_hg38.csv \
  --outdir results/multiple_gwas_eqtl_coloc
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

1. **CLUMP_GWAS**: Process GWAS summary statistics and perform LD-based clumping to identify independent signals
2. **RUN_COLOC**: Perform colocalization analysis between GWAS and eQTL datasets using the coloc R package
3. **COMBINE_COLOC**: Combine all colocalization results and filter based on posterior probability threshold

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--gwas_file` | Single GWAS summary statistics file | null |
| `--gwas_dir` | Directory containing GWAS files | null |
| `--gwas_pattern` | Pattern to match GWAS files | "*.txt" |
| `--gwas_name` | Name for GWAS (only used with single file) | null |
| `--eqtl_file` | Single eQTL data file | null |
| `--eqtl_dir` | Directory containing eQTL files | null |
| `--eqtl_pattern` | Pattern to match eQTL files | "*_cis_MatrixEQTLout.rds" |
| `--gene_location_file` | File with gene coordinates | Required |
| `--pvalue` | P-value threshold for GWAS clumping | 5e-6 |
| `--window_size` | Window size for clumping in base pairs | 1e6 |
| `--coloc_pph4_threshold` | Threshold for PP.H4 to filter results | 0.2 |
| `--outdir` | Output directory | Required |

## Output

The pipeline outputs include:
1. Clumped GWAS data files
2. Colocalization results for each GWAS-eQTL combination
3. Combined and filtered colocalization results

## Citation

If you use this pipeline in your work, please cite:

```
Johnson, Lab et al. (2024). CAUSAL-flow: A Nextflow pipeline for colocalization analysis to identify shared genetic variants. 
GitHub: https://github.com/johnsonlab-ic/CAUSAL-flow
```

## Acknowledgments

This pipeline implements the colocalization method from:

- Giambartolomei C, et al. (2014). Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet. 10(5):e1004383.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

