process run_MR {
    tag "Running MR analysis for ${eqtl_name} and ${gwas_name}"
    label 'process_high'
    
    publishDir "${params.outdir}/mr_results", mode: 'copy'
    
    input:
    path source_R
    path gwas_data
    val gwas_name
    tuple path(eqtl_data), val(eqtl_name)
    val fdr_threshold
    path alleles
    path coloc_results
    
    output:
    path "mr_results_${eqtl_name}_${gwas_name}.txt", emit: results_table
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Source the MR functions
    source("${source_R}")

    # Load required libraries
    library(data.table)
    library(MendelianRandomization)
    library(dplyr)
    
    # Check if this coloc file matches the current GWAS and eQTL
    coloc_filename <- basename("${coloc_results}")
    expected_pattern <- paste0("coloc_results_${eqtl_name}_${gwas_name}")
    
    if (!grepl(expected_pattern, coloc_filename)) {
        cat("This coloc file does not match the current GWAS and eQTL combination. Skipping.\\n")
        quit("no", status = 0)
    }
    
    #read in data
    gwas_data <- readRDS("${gwas_data}")
    eqtl_data <- readRDS("${eqtl_data}")
    allele_df <- data.table::fread("${alleles}")
    
    # Read coloc results and filter for high-confidence colocalization (PP.H4 > 0.8)
    coloc_df <- data.table::fread("${coloc_results}")
    coloc_genes <- coloc_df[coloc_df\$PP.H4 > 0.8,]\$gene
    
    cat(paste0("Found ", length(coloc_genes), " genes with PP.H4 > 0.8 from colocalization results\\n"))
    
    if (length(coloc_genes) == 0) {
        cat("No genes with PP.H4 > 0.8 found. Aborting MR analysis.\\n")
        # Write empty file and exit
        file.create("mr_results_${eqtl_name}_${gwas_name}.txt")
        quit("no", 0)
    }
    
    # Filter eQTL data for coloc genes with significant eQTLs and keep the most significant SNP per gene
    significant_eqtls <- eqtl_data %>%
        filter(gene %in% coloc_genes & FDR < ${fdr_threshold}) %>%
        group_by(gene) %>%
        arrange(FDR) %>%
        slice(1) %>%
        ungroup()
        
    if (nrow(significant_eqtls) == 0) {
        cat("No genes passed both coloc (PP.H4 > 0.8) and eQTL significance (FDR < ${fdr_threshold}) filters. Aborting.\\n")
        file.create("mr_results_${eqtl_name}_${gwas_name}.txt")
        quit("no", 0)
    }
    
    cat(paste0("Running MR on ", nrow(significant_eqtls), " genes that passed both coloc and eQTL filters\\n"))
    
    # Print summary of filtered data
    cat(paste0("Found ", nrow(significant_eqtls), " genes with significant eQTLs (FDR < ${fdr_threshold})\\n"))
    
    # Now, loop through each gene and run MR
    results_list <- list()
    
    # Create output file path
    output_file <- paste0("mr_results_", "${eqtl_name}", "_", "${gwas_name}", ".txt")
    
    # Loop through each gene
    cat("Running MR analysis for each gene:\\n")
    for (i in 1:nrow(significant_eqtls)) {
        gene <- significant_eqtls\$gene[i]
        lead_snp <- significant_eqtls\$SNP[i]
        
        cat(paste0("  Processing gene ", i, "/", nrow(significant_eqtls), ": ", gene, " (SNP: ", lead_snp, ")\\n"))
        
        mr_result <- run_MR_single(
            gwas_data = gwas_data,
            eqtl_data = eqtl_data,
            target_gene = gene,
            allele_df = allele_df,
            lead_snp = lead_snp,
            eqtl_name = "${eqtl_name}", # Add eQTL name parameter
            FDR_cutoff = ${fdr_threshold},
            r2_cutoff = 0.01,
            plink_bin = "/usr/local/bin/plink",
            path_to_binaries = "/app/CAUSAL-flow/EUR"
        )
        
        if (!is.null(mr_result)) {
            results_list[[i]] <- mr_result
        }
    }
    
    # Combine all results
    if (length(results_list) > 0) {
        all_results <- bind_rows(results_list)
        all_results\$gwas_name <- "${gwas_name}"
        all_results\$coloc_status <- "PP.H4_0.8+"
        
        # Write results to file
        write.table(all_results, file = output_file, 
                    sep = "\\t", row.names = FALSE, quote = FALSE)
        
        cat(paste0("MR results saved to ", output_file, "\\n"))
    } else {
        # Write empty results file if no results
        cat("No MR results could be calculated.\\n")
        file.create(output_file)
    }
    """
}

process combine_MR {
    tag "Combining MR results"
    label 'process_low'
    
    publishDir "${params.outdir}/mr_results", mode: 'copy'
    
    input:
    path mr_files
    val pval_threshold
    
    output:
    path "mr_results_combined.txt", emit: combined_results
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    library(data.table)
    
    # List all MR result files
    file_list <- strsplit("${mr_files}", " ")[[1]]
    
    # Read and combine all files
    all_results <- lapply(file_list, function(file) {
      results <- fread(file)
      return(results)
    })
    
    # Combine into one data.frame
    combined_results <- do.call(rbind, all_results)
    
    # Apply threshold filter for significant findings
    filtered_results <- combined_results[combined_results\$IVW_pval < ${pval_threshold}]
    
    # Sort by p-value
    sorted_results <- filtered_results[order(filtered_results\$IVW_pval)]
    
    # Write to file
    write.table(sorted_results, file = "mr_results_combined.txt", 
                sep = "\\t", row.names = FALSE, quote = FALSE)
    """
}
