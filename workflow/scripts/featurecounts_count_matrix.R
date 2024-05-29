#!/usr/bin/Rscript
source("workflow/scripts/script_functions.R", local = TRUE)

if (length(snakemake@input[[1]]) > 1) {
  ids <- read.table(
    snakemake@input[[1]][[1]],
    sep = "\t",
    header = TRUE
  )[, 1]
  ids <- remove_ens_gene_version(ids)

  count_data <- purrr::map_dfc(snakemake@input[[1]], function(x) {
    read.table(
      x,
      sep = "\t",
      header = TRUE
    )[, 7]

    count_data <- cbind.data.frame(ids, count_data)
  })
} else {
  count_data <- read.table(
    snakemake@input[[1]],
    sep = "\t",
    header = TRUE
  )

  count_data <- count_data[, c(1, 7:length(names(count_data)))]
  count_data$Geneid <- remove_ens_gene_version(count_data$Geneid)
}

write.table(
  count_data,
  snakemake@output[[1]],
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
