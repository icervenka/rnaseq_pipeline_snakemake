suppressMessages(library(readr))
suppressMessages(library(purrr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

validate_datasets <- function(dataset_list) {
  # check equal amount of rows
  num_rows <- map(dataset_list, ~nrow)
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
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    dplyr::arrange(`Gene ID`)
})

# check if datasets are valid
validate_datasets(datasets)

combined <- purrr::map_dfr(names(datasets), function(x) {
  datasets[[x]] <- datasets[[x]] %>% mutate(sample = "x")
})

tpm <- combined %>%
  select(sample, `Gene ID`, TPM) %>%
  tidyr::pivot_wider(names_from = "sample", values_from = "TPM")

fpkm <- combined %>%
  select(sample, `Gene ID`, FPKM) %>%
  tidyr::pivot_wider(names_from = "sample", values_from = "FPKM")

# save files
write.table(
  tpm,
  snakemake@output["tpm"],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  fpkm,
  snakemake@output["fpkm"],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  combined,
  snakemake@output["samples_combined"],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
