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

include { PREPARE_GWAS } from './modules/prepare_gwas'
include { PROCESS_EQTL } from './modules/process_eqtl'
include { RUN_COLOC } from './modules/run_coloc'
include { PLOT_RESULTS } from './modules/plot_results'
include { GENERATE_REPORT } from './modules/generate_report'

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Check required parameters
if (!params.gwas_file) {
    log.error "No GWAS file provided. Please specify using --gwas_file"
    exit 1
}

if (!params.eqtl_file) {
    log.error "No eQTL file provided. Please specify using --eqtl_file"
    exit 1
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    // Input files
    gwas_ch = Channel.fromPath(params.gwas_file, checkIfExists: true)
    eqtl_ch = Channel.fromPath(params.eqtl_file, checkIfExists: true)

    // Process GWAS data
    PREPARE_GWAS(gwas_ch)
    prepared_gwas_ch = PREPARE_GWAS.out.gwas_data

    // Process eQTL data
    PROCESS_EQTL(eqtl_ch)
    prepared_eqtl_ch = PROCESS_EQTL.out.eqtl_data

    // Run colocalization analysis
    RUN_COLOC(prepared_gwas_ch, prepared_eqtl_ch)
    coloc_results_ch = RUN_COLOC.out.coloc_results

    // Generate plots
    PLOT_RESULTS(coloc_results_ch)
    plots_ch = PLOT_RESULTS.out.plots

    // Generate final report
    GENERATE_REPORT(coloc_results_ch, plots_ch)
}

/*
========================================================================================
    HELPER FUNCTIONS
========================================================================================
*/

def helpMessage() {
    log.info"""
    COLOC-flow: A Nextflow pipeline for colocalization analysis
    
    Usage:
    nextflow run main.nf --gwas_file <gwas_file> --eqtl_file <eqtl_file>
    
    Required arguments:
      --gwas_file          Path to GWAS summary statistics file
      --eqtl_file          Path to eQTL data file
    
    Optional arguments:
      --output_dir         Directory to save results [default: "${params.output_dir}"]
      --covariates_to_include  Which covariates to include [default: "${params.covariates_to_include}"]
      --help               Show this help message and exit
    """.stripIndent()
}
