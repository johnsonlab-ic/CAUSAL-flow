#!/usr/bin/env nextflow

/*
========================================================================================
    CAUSAL-flow
========================================================================================
    A Nextflow pipeline for colocalization analysis
    Github: https://github.com/johnsonlab-ic/CAUSAL-flow
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
include { run_MR; combine_MR } from './modules/run_MR'
// include { PLOT_RESULTS } from './modules/plot_results'
// include { GENERATE_REPORT } from './modules/generate_report'

/*
========================================================================================
   parameters
========================================================================================
*/
params.outdir="/var/lib/docker/alex_tmp/data/CAUSAL-flow_test/outs"

// GWAS parameters
params.gwas_file="/home/ah3918/GWASi-flow/results/processed/Epilepsy_1_processedGRCh38.txt" // Single GWAS file
params.gwas_dir=null                                                           // Directory with GWAS files (optional)
params.gwas_pattern="*.txt"                                                    // Pattern for GWAS files
params.pvalue=5e-8
params.window_size=1e6
params.coloc_pph4_threshold=0.2 // Threshold for PP.H4 to filter coloc results
params.gwas_name="Epilepsy_1"   // Only used if a single GWAS file is provided

// eQTL parameters
params.eqtl_dir="/var/lib/docker/alex_tmp/data/CAUSAL-flow_test/sc_eqtls/eQTL_outputs" // Directory containing eQTL files
params.eqtl_file=null           // Single eQTL file (optional, takes precedence over dir)
params.eqtl_pattern="*_cis_MatrixEQTLout.rds" // Pattern to match eQTL files
params.gene_location_file="/var/lib/docker/alex_tmp/data/CAUSAL-flow_test/sc_eqtls/expression_matrices/gene_locations.csv"
params.allele_file=null// File with allele information

//MR parameters
params.eqtl_fdr_threshold=0.05 // FDR threshold for eQTLs 
params.mr_pval_threshold=0.05 // P-value threshold for filtering MR results

// add log of input files 
log.info """\
         CAUSAL-flow - Colocalization and Mendelian Randomization Pipeline
         ===================================
         GWAS file    : ${params.gwas_file}
         GWAS dir     : ${params.gwas_dir}
         eQTL file    : ${params.eqtl_file}
         eQTL dir     : ${params.eqtl_dir}
         Output dir   : ${params.outdir}
         eQTL FDR threshold: ${params.eqtl_fdr_threshold}
         MR p-value threshold: ${params.mr_pval_threshold}
         Allele file  : ${params.allele_file}
         """
         .stripIndent()

workflow{
    // Create GWAS channel - either from a single file or a directory
    if (params.gwas_dir) {
        // Multiple GWAS files from directory
        gwas_ch = Channel.fromPath("${params.gwas_dir}/${params.gwas_pattern}")
            .map { file -> 
                def gwas_name = file.baseName
                return tuple(file, gwas_name)
            }
    } else {
        // Single GWAS file (default)
        gwas_ch = Channel.fromPath(params.gwas_file)
            .map { file ->
                return tuple(file, params.gwas_name)
            }
    }
    
    // Create eQTL channel - either from a single file or a directory
    if (params.eqtl_file) {
        // Single eQTL file
        eqtl_ch = Channel.fromPath(params.eqtl_file)
            .map { file -> 
                def eqtl_name = file.name.toString().replaceFirst(/(.+)_cis_MatrixEQTLout\.rds/, '$1')
                return tuple(file, eqtl_name)
            }
    } else {
        // Multiple eQTL files from directory (default)
        eqtl_ch = Channel.fromPath("${params.eqtl_dir}/${params.eqtl_pattern}")
            .map { file -> 
                def eqtl_name = file.name.toString().replaceFirst(/(.+)_cis_MatrixEQTLout\.rds/, '$1')
                return tuple(file, eqtl_name)
            }
    }
    
    // Step 1: Process each GWAS file with clumping
    clumped_gwas_ch = clump_gwas(
        "${baseDir}/R/clump_gwas_source.R",
        gwas_ch.map{it[0]},  // Extract just the file
        gwas_ch.map{params.pvalue},
        gwas_ch.map{params.window_size},
        gwas_ch.map{it[1]}   // Extract just the gwas_name
    )
    
    // Create a combined channel with clumped data and original gwas name
    gwas_data_ch = clumped_gwas_ch.clumped_gwas
        .flatten()
        .map { file -> 
            def gwas_name = file.simpleName.toString().replaceFirst(/(.+)_clumped_gwas/, '$1')
            return tuple(file, gwas_name)
        }
    
    // Step 2: Create combinations of GWAS and eQTL data
    coloc_inputs_ch = gwas_data_ch.combine(eqtl_ch)
    
    // Step 3: Run colocalization analysis for each combination
    coloc_results_ch = run_coloc(
        "${baseDir}/R/run_coloc_functions.R",
        coloc_inputs_ch.map{it[0]},  // gwas_data
        coloc_inputs_ch.map{it[1]},  // gwas_name
        coloc_inputs_ch.map{tuple(it[2], it[3])},  // eqtl tuple
        file(params.gene_location_file)
    )
    
    // Step 4: Combine all colocalization results and filter by PP.H4 threshold
    combine_coloc_ch = combine_coloc(
        coloc_results_ch.results_table.collect(),
        params.coloc_pph4_threshold
    )
    
    // Check if allele file is provided for MR analysis
    if (params.allele_file) {
        log.info "Running Mendelian Randomization analysis..."
        
        // Step 5: Run MR directly using the coloc_complete channel which contains all necessary data
        mr_results_ch = run_MR(
            "${baseDir}/R/run_MR_functions.R",
            coloc_results_ch.coloc_complete.map{it[0]},  // gwas_data
            coloc_results_ch.coloc_complete.map{it[1]},  // gwas_name
            coloc_results_ch.coloc_complete.map{tuple(it[2], it[3])},  // eqtl tuple
            params.eqtl_fdr_threshold,
            file(params.allele_file),
            coloc_results_ch.coloc_complete.map{it[4]}   // coloc results file
        )
        
        // Step 7: Combine all MR results
        combine_MR(
            mr_results_ch.results_table.collect(),
            params.mr_pval_threshold  // P-value threshold for MR results
        )
    } else {
        log.info "Skipping MR analysis as allele file is not provided. Set params.allele_file to run MR."
    }
}
