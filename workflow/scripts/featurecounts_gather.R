#!/usr/bin/Rscript
source("workflow/scripts/script_functions.R", local = TRUE)

print(snakemake@input)
length(snakemake@input[[1]])

if (length(snakemake@input) > 1) {
  ids <- read.table(
    snakemake@input[[1]],
    sep = "\t",
    header = FALSE
  )[, 1]
  ids <- remove_ens_gene_version(ids)

  count_data <- purrr::map_dfc(snakemake@input, function(x) {
    read.table(
      x,
      sep = "\t",
      header = FALSE
    )[, 7]

  })

  count_data <- cbind.data.frame(ids, count_data)
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
  col.names = FALSE,
  row.names = FALSE
)
