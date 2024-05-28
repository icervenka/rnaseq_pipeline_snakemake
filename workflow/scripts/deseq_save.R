#!/usr/bin/Rscript
suppressMessages(library(magrittr))
source("workflow/scripts/script_functions.R", local = TRUE)

# Parse snakemake items --------------------------------------------------------
## input
dds <- readRDS(snakemake@input[["dds"]])
result_array <- readRDS(snakemake@input[["result_array"]])
## output
result_array_ids_out <- snakemake@output[["result_array_ids"]]
sample_expression_out <- snakemake@output[["sample_expression"]]
## params
outdir <- snakemake@params[["outdir"]]
sp_info <- get_species_info(snakemake@params[["species"]])
ids_in <- snakemake@params[["ids_in"]]

# Run --------------------------------------------------------------------------
# create a dataframes with IDs
result_array <- lapply(result_array, function(x) {
  x <- tibble::rownames_to_column(data.frame(x), ids_in)
  x <- add_idcolnames(x, sp_info)
  x <- x[order(x$padj), ]
  return(x)
})

# export deseq corrected sample expression values, merged with ids
sample_expression <- data.frame(DESeq2::counts(dds, normalized = TRUE),
  row.names = rownames(DESeq2::counts(dds))
)

sample_expression <- tibble::rownames_to_column(sample_expression, ids_in) %>%
  add_idcolnames(sp_info)

# fpkm calculation changes the whole data frame, removed until the issue is addressed
# txdb = makeTxDbFromGFF(annotation_gtf, format="gtf")
# txdb_exons = exonsBy(txdb, by="gene")
# txdb_exons = txdb_exons[names(txdb_exons) %in% rownames(dds)]
# rowRanges(dds) = sort(txdb_exons)
# export fpkm values
# fpkm = cbind.data.frame(, fpkm(dds, robust = F))
# fpkm = rownames_to_column(fpkm, 'gene')
# write.table(fpkm, file=snakemake@output[["fpkm"]], quote = F, row.names = F, sep = "\t")

# Save -------------------------------------------------------------------------
write("exporting files", stderr())
# export files with all hits
all <- lapply(names(result_array), function(x) {
  write.table(
    result_array[[x]],
    file = paste0(outdir, x, "_all.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
})

# export significant hits
filtered <- lapply(names(result_array), function(x) {
  write.table(
    result_array[[x]] %>% dplyr::filter(padj < 0.05),
    file = paste0(outdir, x, "_diffexp.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
})

# export normalized expression
write.table(
  sample_expression,
  file = sample_expression_out,
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

saveRDS(result_array, file = result_array_ids_out)
