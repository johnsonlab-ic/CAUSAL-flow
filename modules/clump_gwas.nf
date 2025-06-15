process clump_gwas {
    tag "Clumping ${gwas_name}"
    label 'process_medium'
    
    publishDir "${params.outdir}/clumped_gwas/", mode: 'copy'
    
    input:
    path source_R
    path gwas_file
    val pvalue
    val window_size
    val gwas_name
    
    output:
    path "*clumped_gwas.rds", emit: clumped_gwas
    
    script:

    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(data.table)
    library(ieugwasr)

    work_dir <- getwd()
    temp_dir <- file.path(work_dir, "tmp_dir")
    dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
    Sys.setenv(TMPDIR = temp_dir)
    Sys.setenv(TMP = temp_dir)
    Sys.setenv(TEMP = temp_dir)
    
    #load functions
    source("$source_R")
    # Load example GWAS data
    gwas <- fread("$gwas_file")

    list.files("/app/CAUSAL-flow")
    # Define parameters
    pval <- as.numeric("$pvalue")        
    window <- as.numeric("$window_size")   
    
    # Run region selection
    regions <- select_regions(
        gwas = gwas,
        pval = pval,
        window = window,
        plink_bin = "/usr/local/bin/plink",
        path_to_binaries = "/app/CAUSAL-flow/EUR"
    )
    
    # Save results
    saveRDS(regions, paste0("$gwas_name","_clumped_gwas.rds"))
    
    """
}
