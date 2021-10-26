library(readr)
library(purrr)
library(dplyr)

validate_datasets = function(dataset_list) {
  # check equal amount of rows
  num_rows = map(dataset_list, ~nrow)
  if(length(num_rows %>% unique) != 1) {
    stop("Gene count files from htseq-count don't have equal amount of rows.")
  } else {
    return(TRUE)
  }
}

# read filenames from snakemake and import data
filenames = snakemake@input
datasets = purrr::map(filenames, function(x) {
  readr::read_delim(x,
             '\t',
             col_names = c("geneid", x),
             escape_double = FALSE,
             trim_ws = TRUE) %>%
    dplyr::arrange(geneid)
})

# check if datasets are valid
validate_datasets(datasets)

# extract log information about counting
log_rownames = datasets[[1]] %>%
  dplyr::filter(grepl("__", geneid)) %>%
  dplyr::select(geneid) %>%
  dplyr::transmute(features = gsub("__", "", geneid))
log = purrr::map_dfc(datasets, function(x) {
  x %>%
    dplyr::filter(grepl("__", geneid)) %>%
    dplyr::select(-geneid)
})
log = cbind.data.frame(log_rownames, log)

# extract gene counts
gene_rownames = datasets[[1]] %>%
  dplyr::filter(!grepl("__", geneid)) %>%
  dplyr::select(geneid)
counts = purrr::map_dfc(datasets, function(x) {
  x %>%
    dplyr::filter(!grepl("__", geneid)) %>%
    dplyr::select(-geneid)
})
gene_counts = cbind.data.frame(gene_rownames, counts)

# save files
write.table(gene_counts, snakemake@output['gene_counts'], sep = '\t', quote = F, row.names = F)
write.table(log, snakemake@output['log'], sep = '\t', quote = F, row.names = F)
