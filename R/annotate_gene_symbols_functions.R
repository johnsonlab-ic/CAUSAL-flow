# Adds an Ensembl gene-symbol column to a results table when the "gene" column
# is predominantly Ensembl gene IDs (ENSG*).

library(data.table)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)

#' Add a "symbol" column next to "gene" if >threshold of gene values look like ENSG ids
#'
#' @param input_file Path to a tab-delimited results table with a "gene" column
#' @param output_file Path to write the (possibly annotated) table to
#' @param gene_col Name of the column holding gene identifiers
#' @param threshold Fraction of gene values that must match ^ENSG\\d+ to trigger annotation
annotate_gene_symbols <- function(input_file, output_file, gene_col = "gene", threshold = 0.5) {
  dat <- fread(input_file)

  if (!gene_col %in% names(dat) || nrow(dat) == 0) {
    message("No '", gene_col, "' column or no rows found; writing file through unchanged.")
    fwrite(dat, output_file, sep = "\t")
    return(invisible(dat))
  }

  gene_vals <- dat[[gene_col]]
  gene_vals[is.na(gene_vals)] <- ""
  ensg_frac <- mean(grepl("^ENSG\\d+", gene_vals))

  message(sprintf("%.1f%% of '%s' values look like Ensembl gene IDs.", ensg_frac * 100, gene_col))

  if (ensg_frac <= threshold) {
    message("At or below threshold (", threshold, "); writing file through unchanged.")
    fwrite(dat, output_file, sep = "\t")
    return(invisible(dat))
  }

  message("Above threshold; mapping Ensembl IDs to gene symbols via EnsDb.Hsapiens.v86.")
  ensembl_id <- sub("\\..*$", "", gene_vals)
  symbol <- suppressMessages(mapIds(
    EnsDb.Hsapiens.v86,
    keys = ensembl_id,
    column = "SYMBOL",
    keytype = "GENEID",
    multiVals = "first"
  ))
  symbol <- unname(symbol)
  symbol[is.na(symbol)] <- dat[[gene_col]][is.na(symbol)]

  gene_pos <- which(names(dat) == gene_col)
  dat[, symbol := symbol]
  # names(dat) is re-read post-assignment (rather than reusing a pre-assignment
  # snapshot) because data.table's `:=` mutates column-name attributes by
  # reference, which can silently alias an earlier `names(dat)` capture.
  new_order <- append(setdiff(names(dat), "symbol"), "symbol", after = gene_pos)
  setcolorder(dat, new_order)

  message(sprintf("Resolved %d/%d symbols.", sum(symbol != dat[[gene_col]]), nrow(dat)))

  fwrite(dat, output_file, sep = "\t")
  invisible(dat)
}
