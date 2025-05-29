process run_coloc {
    tag "Running colocalization analysis using cellCOLOC"
    label 'process_high'
    
    publishDir "${params.output_dir}/coloc_results", mode: 'copy'
    
    input:
    path gwas_data
    path eqtl_data
    path gene_loc_file
    path maf_file
    path snp_locs_file
    path indiv_numbers_file
    
    output:
    path "*_cellCOLOC_obj.rds", emit: coloc_obj
    path "*_COLOC_results.txt", emit: coloc_results
    path "*.pdf", optional: true, emit: plots
        
    """
    #!/usr/bin/env Rscript
    
    # Source the cellCOLOC functions
    source("${baseDir}/source_files/cellCOLOC_source.r")
    
    # Load processed data
    gwas <- readRDS("${gwas_data}")
    mateqtlouts <- readRDS("${eqtl_data}")
    
    # Run cellCOLOC analysis
    cellCOLOC(
        GWAS = gwas,
        mateqtlouts = mateqtlouts,
        gene_loc_file = "${gene_loc_file}",
        maf_file = "${maf_file}",
        snp_locs_file = "${snp_locs_file}",
        indiv_numbers_file = "${indiv_numbers_file}",
        preprocess_mateqtlouts = TRUE,
        preprocess_GWAS = TRUE,
        GWASsignif = 5e-8,
        GWAS_window = 1e6,
        GWAS_type = "${gwas_type}",
        GWAS_name = "${gwas_name}",
        cellMR = ${run_cellm
        IVPCA = ${run_ivpca},
        pph4_cutoff = ${pph4_cutoff},
        use_coloc_lead_snp = TRUE
    )
    
    # Generate regional plots if we have colocalized regions
    coloc_results <- read.table("${gwas_name}_COLOC_results.txt", header=TRUE)
    
    # Filter for significant colocalizations
    sig_coloc <- coloc_results[coloc_results\$PP.H4 > ${pph4_cutoff},]
    
    if (nrow(sig_coloc) > 0) {
        # Load the cellCOLOC object
        cellCOLOC_obj <- readRDS("${gwas_name}_cellCOLOC_obj.rds")
        
        # For each significant colocalization, create a regional plot
        for (i in 1:nrow(sig_coloc)) {
            gene <- sig_coloc\$gene[i]
            celltype <- sig_coloc\$cellType[i]
            
            tryCatch({
                plot <- create_regional_association_plot(
                    cellCOLOC_obj,
                    gene = gene,
                    celltype = celltype,
                    gwas_name = "${gwas_name}"
                )
                
                # Save the plot
                pdf(paste0("${gwas_name}_", celltype, "_", gene, "_regional_plot.pdf"), 
                    width = 12, height = 9)
                print(plot)
                dev.off()
            }, error = function(e) {
                message("Error generating plot for gene ", gene, " in celltype ", celltype)
                message(e)
            })
        }
    }
    """
}

