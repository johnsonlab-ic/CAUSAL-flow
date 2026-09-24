process cache_ld_chunks {
    tag "Caching LD chunks for ${gwas_name}"
    label 'process_medium'

    publishDir "${params.outdir}/report_cache/ld_chunks", mode: 'copy'

    input:
    tuple path(gwas_data), val(gwas_name)

    output:
    path "report_cache_ld_${gwas_name}.rds", emit: ld_chunks

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(ieugwasr)

    plink_bin <- "/usr/local/bin/plink"
    ld_bfile <- "/app/CAUSAL-flow/EUR"

    processed_gwas <- readRDS("${gwas_data}")
    if (grepl("chr", processed_gwas\$chr[1]) == FALSE) {
      processed_gwas\$chr <- paste0("chr", processed_gwas\$chr)
    }
    processed_gwas <- split(processed_gwas, processed_gwas\$signif_snp_region)

    # LD depends only on genotype data + which SNPs fall in a region -- never
    # on celltype or gene. Computing it once per region here (in parallel with
    # run_coloc, since this only needs clump_gwas's output) means report
    # generation never has to recompute the same PLINK LD matrix once per
    # celltype for regions that are shared across many celltypes/genes.
    ld_chunks <- list()
    for (region in names(processed_gwas)) {
      snps <- unique(processed_gwas[[region]]\$rsid)
      if (length(snps) < 3) next
      ld <- tryCatch(
        suppressMessages(suppressWarnings(
          ieugwasr::ld_matrix_local(snps, bfile = ld_bfile, plink_bin = plink_bin, with_alleles = FALSE)
        )),
        error = function(e) { message("LD failed for region ", region, ": ", e\$message); NULL }
      )
      if (!is.null(ld)) ld_chunks[[region]] <- ld
    }

    saveRDS(ld_chunks, "report_cache_ld_${gwas_name}.rds")
    message("Cached LD for ", length(ld_chunks), "/", length(processed_gwas), " regions.")
    """
}
