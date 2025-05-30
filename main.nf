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
// include { PROCESS_EQTL } from './modules/process_eqtl'
// include { RUN_COLOC } from './modules/run_coloc'
// include { PLOT_RESULTS } from './modules/plot_results'
// include { GENERATE_REPORT } from './modules/generate_report'

/*
========================================================================================
   parameters
========================================================================================
*/
params.outdir="/var/lib/docker/alex_tmp/data/COLOC-flow_test/outs"

params.gwas_file="/var/lib/docker/alex_tmp/data/COLOC-flow_test/sc_eqtls/GWAS/Epilepsy_1_processedGRCh38.txt"
params.pvalue=5e-8
params.window_size=1e6

workflow{

    clump_gwas(source_R="${baseDir}/R/clump_gwas_source.R",
               gwas_file= "${params.gwas_file}",
               pvalue=params.pvalue,
               window_size=params.window_size)

}
