library(magrittr)
source("workflow/scripts/script_functions.R", local = TRUE)

# Parse snakemake items --------------------------------------------------------
## input
deg_file <- snakemake@input[[1]]
## output

## params
outdir <- snakemake@params[["outdir"]]
fdr <- snakemake@params[["fdr"]]
drop_no_symbol <- snakemake@params[["drop_no_symbol"]]
split_multi_symbols <- snakemake@params[["split_multi_symbols"]]
remove_duplicate_genes <- snakemake@params[["remove_duplicate_genes"]]
# Run --------------------------------------------------------------------------
deg <- read.table(
  deg_file,
  header = TRUE,
  sep = "\t"
) %>%
  dplyr::rename(log2foldchange = "log2.fold_change.")


# move cufflinks gene_ids into symbol column if they don't correspond to known
# genes, if drop_no_symbol is set to true in config, they are removed from the
# data instead
if (drop_no_symbol) {
  deg <- deg %>%
    dplyr::mutate(SYMBOL = gene) %>%
    dplyr::filter(SYMBOL != "-")
} else {
  deg <- deg %>%
    dplyr::mutate(SYMBOL = ifelse(
      gene == "-",
      gene_id,
      gene
    ))
}

# split into separate rows if multiple gene symbols per entry
if (split_multi_symbols) {
  deg <- deg %>%
    tidyr::separate_rows(SYMBOL, sep = ",")
}

# remove duplicate gene names, keeping the ones with lowest adjust p-value
if (remove_duplicate_genes) {
  deg <- deg %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::arrange(q_value, .by_group = TRUE) %>%
    dplyr::slice_head() %>%
    dplyr::arrange(q_value, -abs(log2foldchange)) %>%
    dplyr::ungroup()
} else {
  deg <- deg %>%
    dplyr::arrange(q_value, -abs(log2foldchange))
}

split_deg <- deg %>%
  dplyr::group_by(sample_1, sample_2) %>%
  dplyr::arrange(q_value) %>%
  dplyr::select(
    SYMBOL,
    value_2,
    value_1,
    log2foldchange,
    p_value,
    q_value
  ) %>%
  dplyr::group_map(function(.x, .y) {
    .x <- .x %>%
      dplyr::rename("{.y$sample_1}_value" := value_1) %>% # nolint
      dplyr::rename("{.y$sample_2}_value" := value_2) %>% # nolint
      dplyr::mutate(contrast = paste0(.y$sample_2, "_", .y$sample_1))
  })

# Save -------------------------------------------------------------------------
# files with differential expression
purrr::walk(split_deg, function(x) {
  contrast <- unique(x$contrast)

  x %>%
    dplyr::select(-contrast) %>%
    write.table(
      paste0(outdir, contrast, "_all.txt"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

  x %>%
    dplyr::select(-contrast) %>%
    dplyr::filter(q_value < fdr) %>%
    write.table(
      paste0(outdir, contrast, "_diffexp.txt"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
})

# files with contrasts. Required by the snakemake rule!
deg %>%
  dplyr::select(sample_1, sample_2) %>%
  dplyr::mutate(contrast = paste0(sample_2, "_", sample_1)) %>%
  unique() %>%
  write.table(
    snakemake@output[[1]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )