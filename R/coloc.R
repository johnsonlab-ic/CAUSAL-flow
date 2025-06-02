run_COLOC=function(eqtl_data,gwas_data,gwas_type){


    eqtl_trait = list(
        beta = regionData@beta[[geneName]],  # All values for the specified gene
        varbeta = regionData@std_error[[geneName]]^2,  # Squaring all values for the specified gene
        snp = regionData@beta$SNP,
        type = "quant",
        N = cellCOLOC_obj@cellTypes[[cellTypeName]]@n_indivs,
        MAF = regionData@maf_info$maf,  # Adjusted MAF column name
        sdY=1
    )

}



