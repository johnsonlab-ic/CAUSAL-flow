# Function to select genomic regions from GWAS data
# Expects GWAS data with columns: SNP, CHR, BP, P, etc.

new_ld_clump_local <- function(dat, clump_kb, clump_r2, clump_p, bfile, plink_bin)
{
    # Create a local temp directory
    work_dir <- getwd()
    temp_dir <- file.path(work_dir, "tmp_dir")
    
    # Create the directory if it doesn't exist
    if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
    }
    
    # Make textfile in our custom temp directory
    shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
    random_name <- paste0("clump_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample(1000:9999, 1))
    fn <- file.path(temp_dir, random_name)
    
    # Write the data
    write.table(data.frame(SNP=dat[["rsid"]], P=dat[["pval"]]), 
                file=fn, row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    # Build and run the plink command
    fun2 <- paste0(
        shQuote(plink_bin, type=shell),
        " --bfile ", shQuote(bfile, type=shell),
        " --clump ", shQuote(fn, type=shell), 
        " --clump-p1 ", clump_p, 
        " --clump-r2 ", clump_r2, 
        " --clump-kb ", clump_kb, 
        " --out ", shQuote(fn, type=shell)
    )
    system(fun2)
    
    # Read results and clean up
    clumped_file <- paste(fn, ".clumped", sep="")
    if (file.exists(clumped_file)) {
        res <- tryCatch({
            read.table(clumped_file, header=TRUE)
        }, error = function(e) {
            message("Warning: Could not read clumped file: ", e$message)
            return(data.frame(SNP=character(0)))
        })
    } else {
        message("Warning: Clumped file was not created. Check if plink command ran successfully.")
        return(dat[0,]) # Return empty dataframe with same structure
    }
    
    # Clean up temporary files
    unlink(paste(fn, "*", sep=""))
    
    # Filter results
    y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
    if(nrow(y) > 0) {
        message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel")
    }
    
    return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
}
select_regions = function(
  gwas,                # GWAS data frame with SNP, CHR, BP, P columns
  pval = 5e-8,         # P-value threshold for significant SNPs
  window = 1e6,        # Window size in base pairs
  plink_bin = "/usr/local/bin/plink",    # Will be determined if NULL
  path_to_binaries = "/app/CAUSAL-flow/EUR"  # Path to reference data for clumping
) {
  # Prepare data for ieugwasr (which needs specific column names)
  message("Preparing data for analysis...")
  gwas$rsid <- gwas$SNP  # SNP ID
  gwas$pval <- gwas$P    # P-value
  gwas$pos <- gwas$BP    # Base pair position
  gwas$chr <- gwas$CHR   # Chromosome

  message(plink_bin)
  #now remove thes  
  # Step 1: Split GWAS data by chromosome
  message(paste0(Sys.time(), ": Splitting data by chromosome..."))
  gwas_by_chr <- split(gwas, gwas$CHR)
  
  # Step 2: Filter chromosomes with no significant SNPs
  gwas_by_chr <- lapply(gwas_by_chr, function(chr_data) {
    if (sum(chr_data$P < pval) == 0) {
      return(NULL)
    }
    return(chr_data)
  })
  gwas_by_chr <- Filter(Negate(is.null), gwas_by_chr)
  
  if (length(gwas_by_chr) == 0) {
    message("No regions with significant SNPs found!")
    return(NULL)
  }
  
  # Step 3: Perform LD clumping to identify independent signals
  message(paste0(Sys.time(), ": Identifying independent signals through LD clumping..."))
  kb_window <- window / 1000
  
  # Always use local clumping with reference data from the container
  clumped_snps <- suppressMessages(lapply(
    gwas_by_chr,
    new_ld_clump_local,
    clump_kb = kb_window,
    clump_r2 = 0.001,
    clump_p = pval,
    plink_bin = plink_bin,
    bfile = path_to_binaries
  ))
  
  # Remove failed clumping results
  clumped_snps <- Filter(function(x) !is.null(x) && !all(is.na(x)), clumped_snps)
  
  if (length(clumped_snps) == 0) {
    message("No regions could be clumped successfully!")
    return(NULL)
  }
  
  # Keep only chromosomes with successful clumping
  gwas_by_chr <- gwas_by_chr[names(clumped_snps)]
  
  # Step 4: Define regions around index SNPs
  message(paste0(Sys.time(), ": Defining regions around index SNPs..."))
  
  define_regions <- function(index_snps, chr_data, region_size) {
    # Initialize results
    all_regions <- NULL
    
    # For each index SNP, define a region
    for (i in 1:nrow(index_snps)) {
      snp_pos <- as.numeric(index_snps$pos[i])
      region_start <- snp_pos - region_size
      region_end <- snp_pos + region_size
      
      # Get SNPs within the region
      region_snps <- subset(chr_data, pos >= region_start & pos <= region_end)
      
      if (nrow(region_snps) > 0) {
        # Add region metadata
        region_snps$n_snps_region <- nrow(region_snps)
        region_snps$signif_snp_region <- index_snps$rsid[i]
        
        # Combine with other regions
        all_regions <- rbind(all_regions, region_snps)
      }
    }
    
    return(all_regions)
  }
  
  # Process each chromosome
  half_window <- window / 2
  all_regions <- list()
  
  for (i in 1:length(clumped_snps)) {
    chr_name <- names(clumped_snps)[i]
    all_regions[[i]] <- define_regions(
      clumped_snps[[i]], 
      gwas_by_chr[[chr_name]], 
      half_window
    )
  }
  
  # Combine all regions
  result <- do.call(rbind, all_regions)
  
  message(paste0(Sys.time(), ": Found ", length(unique(result$signif_snp_region)), 
                 " independent regions with ", nrow(result), " total SNPs."))

  result=result %>% select(-SNP, -P, -BP, -CHR)

  return(result)
}


# # # Example usage
# # if (FALSE) {  # Set to TRUE to run the example
#   library(tidyverse)
#   library(data.table)
#   library(ieugwasr)
  
#   # Load example GWAS data
#   gwas <- fread("GWAS/Epilepsy_1_processedGRCh38.txt")
  
#   # Define parameters
#   pval <- 5e-8           # Significance threshold
#   window <- 1e6          # Window size (1Mb)
  
#   # Run region selection
#   regions <- select_regions(
#     gwas = gwas,
#     pval = pval,
#     window = window,
#     plink_bin = genetics.binaRies::get_plink_binary(),
#     path_to_binaries = "/app/CAUSAL-flow/EUR"
#   )
  
#   # Save results
#   saveRDS(regions, "clumped_gwas_regions.rds")
# # }