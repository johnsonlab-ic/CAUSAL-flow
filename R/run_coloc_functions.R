# Functions for running colocalization analysis between GWAS and eQTL data
# 
# These functions will:
# 1. Process each GWAS region
# 2. Find overlapping genes in each region
# 3. Run colocalization for each gene in each region
# 4. Compile results into a comprehensive dataframe

library(tidyverse)
library(data.table)
library(coloc)

#' Run colocalization analysis for a single gene within a GWAS region
#'
#' @param gwas_data GWAS data for a specific region
#' @param eqtl_data eQTL data for the same region
#' @param gene_name Name of the gene to analyze
#' @param trait_type Type of GWAS trait ("quant" or "cc")
#' @param celltype Name of the cell type for eQTL data
#' @param region_name Name of the genomic region
#'
#' @return Named vector with colocalization results
run_single_coloc <- function(gwas_data, eqtl_data, gene_name, trait_type, celltype, region_name) {
  # Filter eQTL data for the specified gene
  gene_eqtl_data <- eqtl_data %>% filter(gene == gene_name)
  
  # If no eQTL data for this gene, return NULL
  if(nrow(gene_eqtl_data) == 0) {
    message(paste0("No eQTL data for gene ", gene_name, " in region ", region_name))
    return(NULL)
  }
  
  # Find common SNPs
  common_snps <- intersect(gwas_data$rsid, gene_eqtl_data$SNP)
  
  # If no common SNPs, return NULL
  if(length(common_snps) == 0) {
    message(paste0("No common SNPs for gene ", gene_name, " in region ", region_name))
    return(NULL)
  }
  
  # Subset to common SNPs
  gwas_subset <- gwas_data[gwas_data$rsid %in% common_snps, ]
  eqtl_subset <- gene_eqtl_data[gene_eqtl_data$SNP %in% common_snps, ]
  
  # Prepare GWAS input for coloc
  coloc_gwas_input <- list(
    beta = as.numeric(gwas_subset$BETA),
    varbeta = as.numeric(gwas_subset$SE)^2,
    snp = as.character(gwas_subset$rsid),
    type = trait_type
  )
  
  # Prepare eQTL input for coloc
  coloc_eqtl_input <- eqtl_subset %>% mutate(se = (beta/t.stat))
  coloc_eqtl_input <- list(
    beta = as.numeric(coloc_eqtl_input$beta),
    varbeta = as.numeric(coloc_eqtl_input$se)^2,
    snp = as.character(coloc_eqtl_input$SNP),
    type = "quant",
    sdY = 1
  )
  
  # Run coloc analysis
  coloc_result <- tryCatch({
    coloc::coloc.abf(coloc_eqtl_input, coloc_gwas_input)
  }, error = function(e) {
    message(paste0("Error in coloc analysis for gene ", gene_name, 
                  " in region ", region_name, ": ", e$message))
    return(NULL)
  })
  
  # If coloc failed, return NULL
  if(is.null(coloc_result)) {
    return(NULL)
  }
  
  # Extract summary results
  resvec <- coloc_result$summary
  
  # Get lead SNP info
  lead_snp_data <- coloc_result$results
  lead_snp_data <- lead_snp_data[order(lead_snp_data$SNP.PP.H4, decreasing=TRUE), ]
  
  if(nrow(lead_snp_data) == 0) {
    message(paste0("No lead SNP found for gene ", gene_name, " in region ", region_name))
    return(NULL)
  }
  
  lead_snp <- lead_snp_data$snp[1]
  snph4 <- lead_snp_data$SNP.PP.H4[1]
  
  # Extract eQTL FDR, eQTL p-value, and GWAS p-value for the lead SNP
  eqtl_info <- eqtl_subset %>% filter(SNP == lead_snp)
  
  if(nrow(eqtl_info) == 0) {
    message(paste0("Lead SNP ", lead_snp, " not found in eQTL data"))
    eqtl_fdr <- NA
    eqtl_pval <- NA
  } else {
    eqtl_fdr <- eqtl_info$FDR[1]
    eqtl_pval <- eqtl_info$p.value[1]
  }
  
  gwas_row <- which(gwas_subset$rsid == lead_snp)
  if(length(gwas_row) == 0) {
    message(paste0("Lead SNP ", lead_snp, " not found in GWAS data"))
    gwas_pval <- NA
  } else {
    gwas_pval <- gwas_subset$pval[gwas_row[1]]
  }
  
  # Combine all results
  full_result <- c(
    region_name, 
    celltype, 
    gene_name, 
    lead_snp, 
    snph4, 
    eqtl_pval, 
    eqtl_fdr, 
    gwas_pval, 
    resvec
  )
  
  names(full_result) <- c(
    "region", 
    "celltype", 
    "gene", 
    "lead_snp",
    "SNP.PP.H4", 
    "eQTL_pval", 
    "eQTL_FDR", 
    "GWAS_pval", 
    "nsnps",
    "PP.H0", 
    "PP.H1", 
    "PP.H2", 
    "PP.H3", 
    "PP.H4"
  )
  
  return(full_result)
}

#' Run colocalization for all genes in a single GWAS region
#'
#' @param region_index Index of the region in the processed_gwas list
#' @param processed_gwas List of processed GWAS data frames by region
#' @param eqtl_list List of eQTL data frames by region
#' @param trait_type Type of GWAS trait ("quant" or "cc")
#' @param celltype Name of the cell type for eQTL data
#'
#' @return Data frame of colocalization results for all genes in the region
run_region_coloc <- function(region_index, processed_gwas, eqtl_list, trait_type, celltype) {
  # Get region data
  region_name <- names(processed_gwas)[region_index]
  gwas_data <- processed_gwas[[region_index]]
  eqtl_data <- eqtl_list[[region_index]]
  
  # Check if eQTL data exists for this region
  if(is.null(eqtl_data) || nrow(eqtl_data) == 0) {
    message(paste0("No eQTL data for region ", region_name))
    return(NULL)
  }
  
  # Get list of genes in the region
  gene_list <- unique(eqtl_data$gene)
  
  if(length(gene_list) == 0) {
    message(paste0("No genes found in region ", region_name))
    return(NULL)
  }
  
  message(paste0("Processing region ", region_name, " with ", length(gene_list), " genes"))
  
  # Find common SNPs between GWAS and eQTL data for this region
  common_snps <- intersect(gwas_data$rsid, eqtl_data$SNP)
  
  # If no common SNPs, skip this region
  if(length(common_snps) == 0) {
    message(paste0("No common SNPs found in region ", region_name))
    return(NULL)
  }
  
  # Subset GWAS and eQTL data to common SNPs
  gwas_data <- gwas_data[gwas_data$rsid %in% common_snps, ]
  eqtl_data <- eqtl_data[eqtl_data$SNP %in% common_snps, ]
  
  # Run coloc for each gene
  all_results <- list()
  for(g in seq_along(gene_list)) {
    gene_name <- gene_list[g]
    
    message(paste0("  Processing gene ", g, "/", length(gene_list), ": ", gene_name))
    
    result <- run_single_coloc(
      gwas_data = gwas_data,
      eqtl_data = eqtl_data,
      gene_name = gene_name,
      trait_type = trait_type,
      celltype = celltype,
      region_name = region_name
    )
    
    if(!is.null(result)) {
      all_results[[gene_name]] <- result
    }
  }
  
  # Combine results into a data frame
  if(length(all_results) > 0) {
    results_df <- as.data.frame(do.call(rbind, all_results))
    # Convert numeric columns to numeric
    numeric_cols <- c("SNP.PP.H4", "eQTL_pval", "eQTL_FDR", "GWAS_pval", 
                     "nsnps", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")
    results_df[numeric_cols] <- lapply(results_df[numeric_cols], as.numeric)
    return(results_df)
  } else {
    return(NULL)
  }
}

#' Run colocalization analysis for all regions and all genes
#'
#' @param processed_gwas List of processed GWAS data frames by region
#' @param eqtl_list List of eQTL data frames by region
#' @param trait_type Type of GWAS trait ("quant" or "cc")
#' @param celltype Name of the cell type for eQTL data
#' @param output_file Optional file path to save the results
#'
#' @return Data frame with all colocalization results
run_all_coloc <- function(processed_gwas, eqtl_list, trait_type, celltype, output_file = NULL) {
  all_region_results <- list()
  
  # Check that processed_gwas and eqtl_list have the same length
  if(length(processed_gwas) != length(eqtl_list)) {
    stop("processed_gwas and eqtl_list must have the same length")
  }
  
  # Process each region
  for(i in seq_along(processed_gwas)) {
    region_results <- run_region_coloc(
      region_index = i,
      processed_gwas = processed_gwas,
      eqtl_list = eqtl_list,
      trait_type = trait_type,
      celltype = celltype
    )
    
    if(!is.null(region_results)) {
      all_region_results[[i]] <- region_results
    }
  }
  
  # Combine results from all regions
  if(length(all_region_results) > 0) {
    final_results <- do.call(rbind, all_region_results)
    
    # Round numeric values for readability
    numeric_cols <- c("SNP.PP.H4", "eQTL_pval", "eQTL_FDR", "GWAS_pval", 
                     "nsnps", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")
    final_results[numeric_cols] <- lapply(final_results[numeric_cols], function(x) {
      if(is.numeric(x)) round(x, 6) else x
    })
    
    # Sort by PP.H4 (posterior probability of colocalization)
    final_results <- final_results[order(final_results$PP.H4, decreasing=TRUE), ]
    
    # Write to file if specified
    if(!is.null(output_file)) {
      write.table(
        final_results, 
        file = output_file, 
        sep = "\t", 
        row.names = FALSE, 
        quote = FALSE
      )
      message(paste0("Results saved to ", output_file))
    }
    
    return(final_results)
  } else {
    message("No colocalization results found")
    return(NULL)
  }
}

get_genes_per_region=function(processed_gwas_region,gene_locs){

  gene_locs<-split(gene_locs,gene_locs$chr)
  gene_locs<-gene_locs[which(names(gene_locs)==unique(processed_gwas_region$chr))]
  gene_locs<-do.call(rbind,gene_locs)
  rownames(gene_locs)<-rep(1:nrow(gene_locs))

  startpos<-min(processed_gwas_region$pos)
  endpos<-max(processed_gwas_region$pos)

  #this is to see whether the end of the gene is within region
  gene_locs_1<-subset(gene_locs,right>startpos & left<endpos)

  #this is to see whether the start of the gene is within region
  gene_locs<-subset(gene_locs,left>startpos & left<endpos)

  #combine the two
  gene_locs_1<-setdiff(gene_locs_1,gene_locs)
  gene_locs<-rbind(gene_locs,gene_locs_1)

  return(gene_locs$geneid)

}

# # Example usage
# if(FALSE) {
#   # Load data
#   gwas_path <- "./GWAS/Epilepsy_1_clumped_gwas.rds"
#   gene_locations_path <- "./expression_matrices/gene_locations.csv"
#   eqtl_path <- "./eQTL_outputs/mOli_cis_MatrixEQTLout.rds"
#   min_snps <- 100
#   trait_type <- "cc"
#   celltype <- "mOli"
#   output_file <- paste0("coloc_results_", celltype, ".txt")
  
#   # Load and process GWAS data
#   processed_gwas <- readRDS(gwas_path) 
#   if(grepl("chr", processed_gwas$chr[1]) == FALSE) {
#     processed_gwas$chr <- paste0("chr", processed_gwas$chr)
#   }
#   processed_gwas <- split(processed_gwas, processed_gwas$signif_snp_region)
  
#   # Load gene locations and find genes in each region
#   gene_locs <- fread(gene_locations_path)
#   genestokeep <- lapply(processed_gwas, get_genes_per_region, gene_locs = gene_locs)
  
#   # Remove regions with no genes
#   regions_no_genes <- which(sapply(genestokeep, length) == 0)
#   if(length(regions_no_genes) > 0) {
#     message(paste0("WARNING: ", length(regions_no_genes), 
#                   " regions have no genes. Removing these regions."))
#     processed_gwas <- processed_gwas[-regions_no_genes]
#     genestokeep <- genestokeep[-regions_no_genes]
#   }
  
#   # Remove small regions
#   processed_gwas <- processed_gwas[sapply(processed_gwas, function(x) {
#     nrow(x) > min_snps
#   })]
  
#   # Load eQTL data and filter by genes in each region
#   eqtls <- readRDS(eqtl_path)
#   eqtl_list <- lapply(genestokeep, function(x) {
#     eqtls[eqtls$gene %in% x, ]
#   })
  
#   # Run colocalization for all regions and all genes
#   coloc_results <- run_all_coloc(
#     processed_gwas = processed_gwas,
#     eqtl_list = eqtl_list,
#     trait_type = trait_type,
#     celltype = celltype,
#     output_file = output_file
#   )
  
#   # Display summary of top hits
#   top_hits <- coloc_results %>% filter(PP.H4 > 0.5)
#   print(paste0("Found ", nrow(top_hits), " colocalization hits with PP.H4 > 0.5"))
#   if(nrow(top_hits) > 0) {
#     print(top_hits)
#   }
# }
