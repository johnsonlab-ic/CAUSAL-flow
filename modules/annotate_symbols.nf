process annotate_gene_symbols {
    tag "Annotating gene symbols for ${result_file.baseName}"
    label 'process_low'
    container 'ghcr.io/johnsonlab-ic/genetics:latest'

    publishDir "${params.outdir}/${publish_subdir}", mode: 'copy', overwrite: true

    input:
    path source_R
    path result_file
    val publish_subdir

    output:
    path "${result_file.name}", emit: annotated_results

    script:
    """
    #!/usr/bin/env Rscript

    source("${source_R}")

    annotate_gene_symbols(
      input_file = "${result_file}",
      output_file = "${result_file.name}"
    )
    """
}
