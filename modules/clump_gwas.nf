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

    custom_tempfile <- function(pattern = "file", tmpdir = NULL, fileext = "") {
    if (is.null(tmpdir)) {
        # Create a local temp directory
        work_dir <- getwd()
        tmpdir <- file.path(work_dir, "tmp_dir")
        
        # Create the directory if it doesn't exist
        if (!dir.exists(tmpdir)) {
        dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
        }
    }
    
    # Generate a random filename in our custom directory
    filename <- paste0(pattern, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", 
                        sample(1000:9999, 1), fileext)
    file.path(tmpdir, filename)
    }
    
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
