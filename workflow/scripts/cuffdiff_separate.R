library(magrittr)
source("workflow/scripts/script_functions.R", local = TRUE)

# Parse snakemake items --------------------------------------------------------
## input
deg_file <- snakemake@input
## output

## params
outdir <- snakemake@params[["outdir"]]
fdr <- snakemake@params[["fdr"]]

# Run --------------------------------------------------------------------------
deg <- read.table(
  deg_file,
  header = TRUE,
  sep = "\t"
)
# Save -------------------------------------------------------------------------

deg %>%
  dplyr::select(sample_1, sample_2) %>%
  dplyr::mutate(contrast = paste0(sample_2, "_", sample_1)) %>%
  write.table(
    snakemake@output[[1]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

deg %>%
  dplyr::group_by(sample_1, sample_2) %>%
  dplyr::arrange(q_value) %>%
  dplyr::select(
    GENE_ID = gene_id,
    SYMBOL = gene,
    value_2,
    value_1,
    log2foldchange = log2.fold_change.,
    p_value,
    q_value
  ) %>%
  dplyr::group_walk(function(.x, .y) {
    .x <- .x %>%
      dplyr::rename("{.y$sample_1}_value" := value_1) %>%
      dplyr::rename("{.y$sample_2}_value" := value_2)

    write.table(
      .x,
      paste0(outdir, .y$sample_2, "_", .y$sample_1, "_all.txt"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }) %>%
  filter(q_value < fdr) %>%
  dplyr::group_walk(function(.x, .y) {
    .x <- .x %>%
      dplyr::rename("{.y$sample_1}_value" := value_1) %>%
      dplyr::rename("{.y$sample_2}_value" := value_2)

    write.table(
      .x,
      paste0(outdir, .y$sample_2, "_", .y$sample_1, "_diffexp.txt"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  })
