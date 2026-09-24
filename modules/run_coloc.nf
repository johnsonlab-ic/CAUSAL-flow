process run_coloc {
    tag "Running colocalization analysis for ${eqtl_name} and ${gwas_name}"
    label 'process_high'
    
    publishDir "${params.outdir}/coloc_results", mode: 'copy', pattern: "coloc_results_*.txt"
    publishDir "${params.outdir}/report_cache/eqtl_chunks", mode: 'copy', pattern: "report_cache_eqtl_*.rds"
    
    input:
    path source_R
    path gwas_data
    val gwas_name
    tuple path(eqtl_data), val(eqtl_name)
    path gene_loc_file
    val do_report
    val eqtl_sdY
    path allele_file, stageAs: 'eqtl_maf.csv'
    path eqtl_summary_file, stageAs: 'eqtl_summary.csv'

    output:
    path "coloc_results_${eqtl_name}_${gwas_name}.txt", emit: results_table
    tuple path(gwas_data), val(gwas_name), path(eqtl_data), val(eqtl_name), path("coloc_results_${eqtl_name}_${gwas_name}.txt"), emit: coloc_complete
    path "report_cache_eqtl_${eqtl_name}_${gwas_name}.rds", optional: true, emit: report_cache

        
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

    # Colocalization needs the full set of tested SNP-gene associations per region
    # (coloc.abf compares H3 vs H4 across all tested SNPs); an FDR/p-value-filtered
    # eQTL file silently starves regions down to a handful of SNPs and produces
    # unreliable PP.H4 estimates. Flag the common case where every row already
    # cleared a significance filter.
    if (nrow(eqtls) > 0 && "FDR" %in% names(eqtls) && all(eqtls\$FDR <= 0.05, na.rm = TRUE)) {
      cat(paste0(
        "WARNING: eQTL input '${eqtl_name}' looks FDR-filtered (every row has FDR <= 0.05). ",
        "coloc requires the FULL (unfiltered) set of tested SNP-gene associations per region, ",
        "not just significant hits, or PP.H4 estimates will be unreliable, especially where nsnps is low. ",
        "Point --eqtl_dir/--eqtl_pattern at the unfiltered eQTL output instead.\\n"
      ))
    }

    eqtl_list <- lapply(genestokeep, function(x) {
      eqtls[eqtls\$gene %in% x, ]
    })

    # Optional: cache the per-region eQTL subset already computed above, for
    # later report generation. Costs nothing extra -- eqtl_list is already in
    # memory here, so this avoids a separate process re-reading the (huge)
    # raw eQTL file just to redo the same region/gene filtering.
    if (toupper("${do_report}") == "TRUE") {
      saveRDS(eqtl_list, "report_cache_eqtl_${eqtl_name}_${gwas_name}.rds")
    }
    
    # eQTL sdY: "estimate" (from MAF + N) or a fixed number. MAF and N are also
    # used for the mismatch guardrail when a fixed number is given. Empty
    # staged files mean the input wasn't supplied (see NO_FILE in main.nf).
    eqtl_sdY_mode <- "${eqtl_sdY}"
    if (eqtl_sdY_mode != "estimate") eqtl_sdY_mode <- as.numeric(eqtl_sdY_mode)

    eqtl_n <- NA_real_
    if (file.size("eqtl_summary.csv") > 0) {
      eqtl_summary <- fread("eqtl_summary.csv")
      eqtl_n <- as.numeric(eqtl_summary[celltype == eqtl_name]\$n_individuals[1])
    }

    eqtl_maf <- NULL
    if (file.size("eqtl_maf.csv") > 0) {
      maf_dt <- fread("eqtl_maf.csv", select = c("snp", "maf"))
      maf_dt <- maf_dt[snp %in% unlist(lapply(processed_gwas, function(x) x\$rsid))]
      eqtl_maf <- setNames(maf_dt\$maf, maf_dt\$snp)
    }

    if (identical(eqtl_sdY_mode, "estimate") && (is.na(eqtl_n) || is.null(eqtl_maf))) {
      stop(paste0("--eqtl_sdY estimate needs N and MAF for celltype '", eqtl_name, "': ",
                  "not found in --eqtl_summary_file (n_individuals) and/or --allele_file (maf)."))
    }
    if (is.na(eqtl_n) || is.null(eqtl_maf)) {
      cat("WARNING: no N and/or MAF for '${eqtl_name}'; sdY mismatch guardrail disabled.\\n")
    }

    sdy_opts <- list(mode = eqtl_sdY_mode, n = eqtl_n, maf = eqtl_maf)

    # Set output files
    output_file <- paste0("coloc_results_", eqtl_name, "_${gwas_name}.txt")


    # Run colocalization
    coloc_results <- run_all_coloc(
      processed_gwas = processed_gwas,
      eqtl_list = eqtl_list,
      trait_type = trait_type,
      celltype = eqtl_name,
      sdy_opts = sdy_opts,
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
    
    # Collect staged coloc result files directly from the task working directory.
    # This avoids breaking paths when file names contain spaces.
    file_list <- list.files(pattern = "^coloc_results_.*[.]txt", full.names = TRUE)
    if (length(file_list) == 0) {
      stop("No coloc result files were staged for combine_coloc")
    }
    
    # Read and combine all files
    all_results <- lapply(file_list, function(file) {
      results <- fread(file)
      # Make sure GWAS column exists
      return(results)
    })
    
    # Combine into one data.frame
    combined_results <- do.call(rbind, all_results)
    
    # Apply threshold filter
    filtered_results <- combined_results
    
    # Sort by PP.H4
    sorted_results <- filtered_results[order(-filtered_results\$PP.H4)]
    
    # Write to file
    write.table(sorted_results, file = "coloc_results_combined.txt", 
                sep = "\\t", row.names = FALSE, quote = FALSE)
    """
}



