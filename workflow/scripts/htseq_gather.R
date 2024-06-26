#!/usr/bin/Rscript
library(magrittr)
source("workflow/scripts/script_functions.R", local = TRUE)

validate_datasets <- function(dataset_list) {
  # check equal amount of rows
  num_rows <- purrr::map(dataset_list, ~nrow)
  if (length(num_rows %>% unique()) != 1) {
    stop("Gene count files from htseq-count don't have equal amount of rows.")
  } else {
    return(TRUE)
  }
}

# read filenames from snakemake and import data
filenames <- snakemake@input
datasets <- purrr::map(filenames, function(x) {
  readr::read_delim(x,
    "\t",
    col_names = c("geneid", x),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    dplyr::arrange(geneid)
})

# check if datasets are valid
validate_datasets(datasets)

if (snakemake@params[["kind"]] == "count") {
  # extract gene counts
  gene_rownames <- datasets[[1]] %>%
    dplyr::filter(!grepl("__", geneid)) %>%
    dplyr::pull(geneid)
  gene_rownames <- remove_ens_gene_version(gene_rownames)

  counts <- purrr::map_dfc(datasets, function(x) {
    x %>%
      dplyr::filter(!grepl("__", geneid)) %>%
      dplyr::select(-geneid)
  })
  gene_counts <- cbind.data.frame(GeneId = gene_rownames, counts)

  # save count file
  write.table(
    gene_counts,
    snakemake@output[[1]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
} else if (snakemake@params[["kind"]] == "log") {
  # extract log information about counting
  log_rownames <- datasets[[1]] %>%
    dplyr::filter(grepl("__", geneid)) %>%
    dplyr::select(geneid)

  log <- purrr::map(datasets, function(x) {
    x %>%
      dplyr::filter(grepl("__", geneid)) %>%
      dplyr::select(-geneid) %>%
      dplyr::mutate(features = log_rownames$geneid) %>%
      dplyr::relocate(features)
  })
  # save log files
  purrr::map(log, function(x) {
    write.table(
      x,
      snakemake@output[[1]],
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  })
}
