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


run_eQTL_GWAS_COLOC=function(processed_gwas,eqtl,gene_locs){

   processed_gwas<-split(processed_gwas,processed_gwas$chr)
   genestokeep=lapply(processed_gwas,get_genes_per_region,gene_locs=gene_locs)






}

#### test code in docker
gwas_path="./GWAS/Epilepsy_1_clumped_gwas.rds"
gene_locations_path="./expression_matrices/gene_locations.csv"
eqtl_path="./eQTL_outputs/mOli_cis_MatrixEQTLout.rds"
min_snps=100
trait_type="cc"
celltype="mOli"


library(tidyverse)
library(data.table)
processed_gwas=readRDS(gwas_path) 
if(grepl("chr",processed_gwas$chr[1])==FALSE){
  processed_gwas$chr<-paste0("chr",processed_gwas$chr)
}
processed_gwas=split(processed_gwas,processed_gwas$signif_snp_region)
gene_locs=fread(gene_locations_path)
genestokeep=lapply(processed_gwas,get_genes_per_region,gene_locs=gene_locs)

# remove regions with no genes.
regions_no_genes=which(sapply(genestokeep,length)==0)
if(length(regions_no_genes)>0){
  message(paste0("WARNING: ",length(regions_no_genes)," regions have no genes. Removing these regions."))
  processed_gwas=processed_gwas[-regions_no_genes]
  genestokeep=genestokeep[-regions_no_genes]
}

#remove small regions, and remove list of empty elements
processed_gwas=processed_gwas[sapply(processed_gwas,function(x){nrow(x)>min_snps})]

# now the eqtls 
eqtls=readRDS(eqtl_path)
eqtl_list=lapply(genestokeep,function(x){
  eqtls[eqtls$gene %in% x,]
})

i=1
regionName=names(processed_gwas)[i]
gwas_data=processed_gwas[[i]]
eqtl_data=eqtl_list[[i]]

# make sure the same SNPs in same order
common_snps=intersect(gwas_data$rsid,eqtl_data$SNP)
# message to see how many were lost and which SNPs were removed
lost_gwas_snps <- setdiff(gwas_data$rsid, common_snps)
lost_eqtl_snps <- setdiff(eqtl_data$SNP, common_snps)
if(length(lost_gwas_snps) > 0 | length(lost_eqtl_snps) > 0){
    message(paste0("WARNING: SNPs were lost in region ", i, "."))
    if(length(lost_gwas_snps) > 0){
        message(paste0("  GWAS SNPs removed: ", length(lost_gwas_snps)))
    }
    if(length(lost_eqtl_snps) > 0){
        message(paste0("  eQTL SNPs removed: ", length(lost_gwas_snps)))
    }
}

gwas_data=gwas_data[gwas_data$rsid %in% common_snps,]
eqtl_data=eqtl_data[eqtl_data$SNP %in% common_snps,]

coloc_gwas_input=list(beta=as.numeric(gwas_data$BETA),
varbeta=as.numeric(gwas_data$SE)^2,
snp=as.character(gwas_data$rsid),
type=trait_type)



###now create COLOC inputs
g=1
genelist=unique(eqtl_data$gene)
geneName=genelist[g]
coloc_eqtl_input=eqtl_data %>% filter(gene==geneName) %>% mutate(se=(beta/t.stat))
coloc_eqtl_input=list(beta=as.numeric(coloc_eqtl_input$beta),
                      varbeta=as.numeric(coloc_eqtl_input$se)^2,
                      snp=as.character(coloc_eqtl_input$SNP),
                      type="quant",
                      sdY=1)

coloc_result <-coloc::coloc.abf(coloc_eqtl_input, coloc_gwas_input)
resvec=coloc_result$summary
lead_snp_data = coloc_result$results
lead_snp_data = lead_snp_data[order(lead_snp_data$SNP.PP.H4, decreasing=T), ]

lead_snp = lead_snp_data$snp[1]
snph4 = lead_snp_data$SNP.PP.H4[1]

# Extract eQTL FDR, eQTL p-value, and GWAS p-value for the lead SNP
eqtl_fdr = eqtl_data %>% filter(gene==geneName,SNP==lead_snp) %>% pull(FDR)
eqtl_pval =  eqtl_data %>% filter(gene==geneName,SNP==lead_snp) %>% pull(p.value)
gwas_pval = gwas_data$pval[gwas_data$rsid == lead_snp]

# Append these values to resvec
resvec = c(regionName, celltype, geneName, lead_snp, snph4, eqtl_pval, eqtl_fdr, gwas_pval, resvec)
names(resvec) = c("region", "celltype", "gene", "lead_snp","SNP.PP.H4", "eQTL_pval", "eQTL_FDR", "GWAS_pval", "nsnps","PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")

