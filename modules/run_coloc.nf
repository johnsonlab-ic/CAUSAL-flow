process run_coloc {
    tag "Running colocalization analysis for ${eqtl_name} and ${gwas_name}"
    label 'process_high'
    
    publishDir "${params.outdir}/coloc_results", mode: 'copy', pattern: "coloc_results_*.txt"
    
    input:
    path source_R
    path gwas_data
    val gwas_name
    tuple path(eqtl_data), val(eqtl_name)
    path gene_loc_file
    
    output:
    path "coloc_results_${eqtl_name}_${gwas_name}.txt", emit: results_table
    tuple path(gwas_data), val(gwas_name), path(eqtl_data), val(eqtl_name), path("coloc_results_${eqtl_name}_${gwas_name}.txt"), emit: coloc_complete

        
    script:
    
    """
    #!/usr/bin/env Rscript
    
    # Source the colocalization functions
    source("${source_R}")
    
    # Set defaults
    trait_type <- "cc"  # Assuming all trait types are case-control 
    min_snps <- 100     # Minimum SNPs per region set to 100
    eqtl_name <- "${eqtl_name}"
    gwas_name <- "${gwas_name}"


    # Load required libraries
    library(tidyverse)
    library(data.table)
    
    # Load and process GWAS data
    processed_gwas <- readRDS("${gwas_data}") 
    if(grepl("chr", processed_gwas\$chr[1]) == FALSE) {
      processed_gwas\$chr <- paste0("chr", processed_gwas\$chr)
    }
    processed_gwas <- split(processed_gwas, processed_gwas\$signif_snp_region)
    
    # Load gene locations and find genes in each region
    gene_locs <- fread("${gene_loc_file}")
    
    # Find genes in each region
    genestokeep <- lapply(processed_gwas, get_genes_per_region, gene_locs = gene_locs)
    
    # Remove regions with no genes
    regions_no_genes <- which(sapply(genestokeep, length) == 0)
    if(length(regions_no_genes) > 0) {
      cat(paste0("WARNING: ", length(regions_no_genes), " regions have no genes. Removing.\\n"))
      processed_gwas <- processed_gwas[-regions_no_genes]
      genestokeep <- genestokeep[-regions_no_genes]
    }
    
    # Remove small regions
    processed_gwas <- processed_gwas[sapply(processed_gwas, function(x) {
      nrow(x) > min_snps
    })]
    
    # Load eQTL data
    eqtls <- readRDS("${eqtl_data}")
    eqtl_list <- lapply(genestokeep, function(x) {
      eqtls[eqtls\$gene %in% x, ]
    })
    
    # Set output files
    output_file <- paste0("coloc_results_", eqtl_name, "_${gwas_name}.txt")
    
    
    # Run colocalization
    coloc_results <- run_all_coloc(
      processed_gwas = processed_gwas,
      eqtl_list = eqtl_list,
      trait_type = trait_type,
      celltype = eqtl_name,
      output_file = output_file
    )
    coloc_results\$GWAS=gwas_name
    # Save results
    write.table(coloc_results, file = output_file, sep = "\\t", row.names = FALSE, quote = FALSE)
    

    """
}

process combine_coloc {
    tag "Combining colocalization results"
    label 'process_low'
    
    publishDir "${params.outdir}/coloc_results", mode: 'copy'
    
    input:
    path coloc_files
    val pp_threshold
    
    output:
    path "coloc_results_combined.txt", emit: combined_results
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    library(data.table)
    
    # List all coloc result files
    file_list <- strsplit("${coloc_files}", " ")[[1]]
    
    # Read and combine all files
    all_results <- lapply(file_list, function(file) {
      results <- fread(file)
      # Make sure GWAS column exists
      return(results)
    })
    
    # Combine into one data.frame
    combined_results <- do.call(rbind, all_results)
    
    # Apply threshold filter
    filtered_results <- combined_results[combined_results\$PP.H4 > ${pp_threshold}]
    
    # Sort by PP.H4
    sorted_results <- filtered_results[order(-filtered_results\$PP.H4)]
    
    # Write to file
    write.table(sorted_results, file = "coloc_results_combined.txt", 
                sep = "\\t", row.names = FALSE, quote = FALSE)
    """
}



