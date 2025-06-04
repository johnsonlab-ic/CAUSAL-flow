#!/usr/bin/env nextflow

/*
========================================================================================
    COLOC-flow
========================================================================================
    A Nextflow pipeline for colocalization analysis
    Github: https://github.com/yourusername/COLOC-flow
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { clump_gwas } from './modules/clump_gwas'
include { run_coloc; combine_coloc } from './modules/run_coloc'
// include { PLOT_RESULTS } from './modules/plot_results'
// include { GENERATE_REPORT } from './modules/generate_report'

/*
========================================================================================
   parameters
========================================================================================
*/
params.outdir="/var/lib/docker/alex_tmp/data/COLOC-flow_test/outs"

// GWAS parameters
params.gwas_file="/home/ah3918/GWASi-flow/results/processed/Epilepsy_1_processedGRCh38.txt"
params.pvalue=5e-8
params.window_size=1e6
params.coloc_pph4_threshold=0.5 // Threshold for PP.H4 to filter coloc results
params.gwas_name="Epilepsy_1"

// eQTL parameters
params.eqtl_dir="/var/lib/docker/alex_tmp/data/COLOC-flow_test/sc_eqtls/eQTL_outputs" // Directory containing eQTL files
params.eqtl_pattern="*_cis_MatrixEQTLout.rds" // Pattern to match eQTL files
params.gene_location_file="/var/lib/docker/alex_tmp/data/COLOC-flow_test/sc_eqtls/expression_matrices/gene_locations.csv"


// add log of input files 
log.info """\
         COLOC-flow - Colocalization Analysis Pipeline
         ===================================
         GWAS file    : ${params.gwas_file}
         eQTL dir     : ${params.eqtl_dir}
         Output dir   : ${params.outdir}
         GWAS name    : ${params.gwas_name}
         """
         .stripIndent()

workflow{
    // Step 1: GWAS clumping
    clump_gwas(
        source_R="${baseDir}/R/clump_gwas_source.R",
        gwas_file=params.gwas_file,
        pvalue=params.pvalue,
        window_size=params.window_size,
        gwas_name=params.gwas_name
    )
    
    // Create a channel for eQTL files with extracted cell/tissue type names
    Channel
        .fromPath("${params.eqtl_dir}/${params.eqtl_pattern}")
        .map { file -> 
            // Extract the eQTL name from the filename (e.g., "ChoroidPlexus" from "ChoroidPlexus_cis_MatrixEQTLout.rds")
            def eqtl_name = file.name.toString().replaceFirst(/(.+)_cis_MatrixEQTLout\.rds/, '$1')
            return tuple(file, eqtl_name)
        }
        .set { eqtl_files_ch }
    
    // Step 2: Run colocalization analysis for each eQTL file
    run_coloc(
        "${baseDir}/R/run_coloc_functions.R",
        clump_gwas.out.clumped_gwas,
        params.gwas_name,
        eqtl_files_ch,
        file(params.gene_location_file)
    )
    
    // Step 3: Combine all colocalization results and filter by PP.H4 threshold
    combine_coloc(
        run_coloc.out.results_table.collect(),
       params.coloc_pph4_threshold // PP.H4 threshold
    )
}
