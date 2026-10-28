# CAUSAL-flow

A Nextflow pipeline for colocalization (`coloc`) and Mendelian Randomization analysis between GWAS and eQTL datasets.

## Overview

This pipeline implements colocalization analysis using the `coloc` R package to identify genetic variants that influence both disease risk and gene expression. The workflow clumps GWAS summary statistics into independent loci, runs coloc for every GWAS × eQTL (cell type) pair, and optionally runs MR on colocalized genes.

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
  --eqtl_file <path/to/<celltype>_cis_MatrixEQTLout.rds> \
  --gene_location_file <path/to/gene_locations.csv> \
  --allele_file <path/to/MAF_mat.csv> \
  --eqtl_summary_file <path/to/eqtl_summary.csv> \
  --outdir <output_directory>
```

By default the eQTL sdY is estimated (`--eqtl_sdY estimate`), which requires `--allele_file` and `--eqtl_summary_file`. For standardised eQTL inputs, pass `--eqtl_sdY 1` instead and both become optional (without `--allele_file`, MR is skipped). The examples below omit these for brevity.

To ensure you always run the latest version, before each run use;

```bash
nextflow drop johnsonlab-ic/CAUSAL-flow
```

### Input Files

1. **GWAS Summary Statistics** (e.g. GWASi-flow output):
   - Columns: `SNP`, `CHR`, `BP`, `P`, `BETA`, `SE`; MR additionally uses `A1`, `A2`
   - Coordinates must be on the same build as the gene location file

2. **eQTL Data** (one file per cell type, e.g. sceQTL-flow / MatrixEQTL output):
   - RDS data frame with columns: `SNP`, `gene`, `beta`, `t.stat`, `p.value`, `FDR`
   - Must be the full, unfiltered cis output (coloc needs all tested SNP-gene pairs)
   - Cell type name is taken from the filename: `<celltype>_cis_*.rds`

3. **Gene Location File**:
   - CSV with columns: `geneid`, `chr` (e.g. `chr1`), `left`, `right`; `geneid` must match the eQTL `gene` column

4. **Allele File** (`--allele_file`, optional): sceQTL-flow `MAF_mat.csv` with columns `snp`, `ref`, `alt`, `maf`

5. **eQTL Summary File** (`--eqtl_summary_file`, optional): sceQTL-flow `eqtl_summary.csv` with columns `celltype`, `n_individuals`

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

- `standard`: Default profile; runs locally with Docker
- `offline`: Runs locally with Docker, 32 GB per task, local executor capped at 300 GB total
- `imperial`: PBS Pro job scheduler with Singularity containers (Imperial HPC)

Run with a specific profile:

```bash
nextflow run main.nf -profile imperial --gwas_file <path> --eqtl_file <path>
```

## Pipeline Steps

1. **CLUMP_GWAS**: Resolve duplicated (multiallelic) rsIDs, then perform LD-based clumping to identify independent signals
2. **RUN_COLOC**: Perform colocalization analysis between GWAS and eQTL datasets using the coloc R package
3. **COMBINE_COLOC**: Combine all colocalization results, sorted by PP.H4, and add gene symbols/biotypes if genes are Ensembl IDs
4. **RUN_MR** / **COMBINE_MR** (only with `--allele_file`): Wald-ratio MR on the coloc lead SNP for genes with PP.H4 > 0.8, eQTL FDR < `--eqtl_fdr_threshold`; combined results filtered at `--mr_pval_threshold`
5. **CACHE_LD_CHUNKS** (only with `--report true`): caches per-region LD matrices and eQTL subsets for report generation

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
| `--pvalue` | P-value threshold for GWAS clumping | 5e-8 |
| `--window_size` | Window size for clumping in base pairs | 1e6 |
| `--coloc_pph4_threshold` | Threshold for PP.H4 to filter results | 0.2 |
| `--allele_file` | In-sample MAF file (sceQTL-flow `MAF_mat.csv`); used for eQTL sdY estimation, MR, and resolving multiallelic GWAS rsIDs (keeps the copy matching eQTL alleles; without it, all duplicated rsIDs are dropped) | null |
| `--eqtl_summary_file` | sceQTL-flow `eqtl_summary.csv`; per-celltype N (`n_individuals`) for eQTL sdY estimation | null |
| `--eqtl_sdY` | eQTL phenotype SD for coloc: `estimate` (from MAF + N, correct at any scale) or a number (e.g. `1` for standardised eQTLs). A fixed value is still checked against the estimate and flagged if they differ >2x | estimate |
| `--eqtl_fdr_threshold` | eQTL FDR threshold for MR instruments | 0.05 |
| `--mr_pval_threshold` | P-value threshold for combined MR results | 0.05 |
| `--report` | Cache LD/eQTL chunks for report generation | false |
| `--outdir` | Output directory | Required |

## Output

The pipeline outputs include:
1. `clumped_gwas/` — clumped GWAS regions (RDS)
2. `coloc_results/` — per GWAS × cell type results and `coloc_results_combined.txt`
3. `mr_results/` — per pair and combined MR results (if `--allele_file`)
4. `report_cache/` — LD and eQTL chunks (if `--report true`)

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

