# Functions for running Mendelian Randomization analysis between GWAS and eQTL data
#
# These functions will:
# 1. Process each GWAS and identify suitable instruments
# 2. Harmonize GWAS and eQTL data
# 3. Run MR analyses with multiple methods
# 4. Compile results into a comprehensive dataframe

library(tidyverse)
library(data.table)
library(MendelianRandomization)
library(ieugwasr)



run_MR_single=function(gwas_data=NULL,
eqtl_data=NULL,
target_gene=NULL,
allele_df=NULL,
lead_snp=NULL,
eqtl_name=NULL,
FDR_cutoff=0.05,
r2_cutoff=0.01,
plink_bin,
path_to_binaries){

    # Filter eQTL data for the target gene and FDR threshold
    eqtl_data=eqtl_data %>% filter(gene==target_gene, FDR<=FDR_cutoff)
    if(nrow(eqtl_data)==0){
        message(paste0(Sys.time(),": No eQTLs found for gene ",target_gene," with FDR <= ",FDR_cutoff))
        return(NULL)
    }
    
    # Filter to only the lead SNP and calculate standard error
    eqtl_data=eqtl_data %>% filter(SNP==lead_snp) %>% mutate(se = abs(beta/t.stat))
    if(nrow(eqtl_data)==0){
        message(paste0(Sys.time(),": Lead SNP ",lead_snp," not found for gene ",target_gene))
        return(NULL)
    }
    
    # Join with allele information
    eqtl_data=eqtl_data %>% left_join(allele_df, by=c("SNP"="snp")) 
    eqtl_data = eqtl_data %>% rename(eqtl_effect_allele=ref, eqtl_other_allele=alt)
    
    # Filter GWAS data for the lead SNP using rsid column
    gwas_data = gwas_data %>% filter(rsid == lead_snp)
    if(nrow(gwas_data) == 0){
        message(paste0(Sys.time(), ": Lead SNP ", lead_snp, " not found in GWAS data"))
        return(NULL)
    }
    
    # Use A1 and A2 as effect alleles
    gwas_data = gwas_data %>% rename(gwas_effect_allele = A2, gwas_other_allele = A1)

    # Create MR dataframe
    mr_df = data.frame(
        SNP = eqtl_data$SNP,
        eqtl_effect_allele = eqtl_data$eqtl_effect_allele,
        eqtl_other_allele = eqtl_data$eqtl_other_allele,
        gwas_effect_allele = gwas_data$gwas_effect_allele,
        gwas_other_allele = gwas_data$gwas_other_allele,
        eqtl_beta = eqtl_data$beta,
        eqtl_se = eqtl_data$se,
        gwas_beta = gwas_data$BETA,
        gwas_se = gwas_data$SE
    )

    # Harmonize alleles (flip beta if allele orientations don't match)
    mr_df=mr_df %>%
    mutate(
        eqtl_beta = ifelse(eqtl_effect_allele != gwas_effect_allele & eqtl_effect_allele == gwas_other_allele, -eqtl_beta, eqtl_beta)
    )

    # Run MR analysis
    tryCatch({
        MRInputObject <- MendelianRandomization::mr_input(
            bx=mr_df$eqtl_beta,
            bxse=mr_df$eqtl_se,
            by=mr_df$gwas_beta,
            byse=mr_df$gwas_se,
            snps=mr_df$SNP
        )
        
        res=MendelianRandomization::mr_ivw(MRInputObject)
        res=data.frame(Method="IVW",Estimate=res@Estimate,`P-value`=res@Pvalue,check.names=F,std_error=res@StdError)
        
        mr_results=res %>% filter(Method == "IVW") %>%
        mutate(
            IVW_beta = signif(Estimate, digits = 4),
            IVW_pval = signif(`P-value`, digits = 4),
            IVW_se=signif(std_error, digits = 4)
        ) %>%
        select(IVW_beta, IVW_pval,IVW_se)
        
        mr_results$IVs <- paste(unique(mr_df$SNP), collapse = ",")
        mr_results$gene <- target_gene  # Use the provided target_gene parameter
        
        # Handle celltype more robustly
        if("celltype" %in% colnames(eqtl_data) && !is.null(eqtl_data$celltype[1])) {
            mr_results$celltype <- eqtl_data$celltype[1]
        } else {
            # Use the eqtl_name parameter if provided, otherwise use a default
            mr_results$celltype <- if(!is.null(eqtl_name)) eqtl_name else "unknown"
        }
        
        mr_results <- mr_results %>% select(celltype, gene, IVs, IVW_beta, IVW_pval, IVW_se)
        return(mr_results)
    }, error = function(e) {
        message(paste0(Sys.time(),": Error in MR analysis for gene ",target_gene,": ", e$message))
        return(NULL)
    })
}