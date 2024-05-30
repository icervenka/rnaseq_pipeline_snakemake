#!/usr/bin/Rscript
library(magrittr)
source("workflow/scripts/script_functions.R", local = TRUE)

# read filenames from snakemake and import data
filenames <- snakemake@input
datasets <- purrr::map(filenames, function(x) {
  readr::read_delim(x,
    "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    dplyr::arrange(`Gene ID`)
}) %>%
  setNames(snakemake@params[["samples"]])

combined <- purrr::map_dfr(names(datasets), function(x) {
  datasets[[x]] %>% mutate(sample = x)
}) %>%
  distinct(`Gene ID`, sample, .keep_all = TRUE)

tpm <- combined %>%
  dplyr::select(sample, `Gene ID`, TPM) %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = "TPM",
    values_fill = NA
  )

fpkm <- combined %>%
  dplyr::select(sample, `Gene ID`, FPKM) %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = "FPKM",
    values_fill = NA
  )

# export data
write.table(
  tpm,
  snakemake@output[["tpm"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  fpkm,
  snakemake@output[["fpkm"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  combined,
  snakemake@output[["samples_combined"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
