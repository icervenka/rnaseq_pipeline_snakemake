#!/usr/bin/Rscript
suppressMessages(library(magrittr))
source("workflow/scripts/script_functions.R", local = TRUE)

result_array <- readRDS(snakemake@input[["result_array"]])
dds <- readRDS(snakemake@input[["dds"]])
outdir <- snakemake@params[["outdir"]]
result_array_ids <- snakemake@output[["result_array_ids"]]
sample_expression <- snakemake@output[["sample_expression"]]
sp_arr <- get_species_info(snakemake@params[["species"]])
ids_in <- snakemake@params[["gene_ids_in"]]

# create a dataframes with IDs
result_array <- lapply(result_array, function(x) {
  x <- tibble::rownames_to_column(data.frame(x), ids_in)
  x <- add_idcolnames(x, sp_arr)
  x <- x[order(x$padj), ]
  return(x)
})

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

# fpkm calculation changes the whole data frame, removed until the issue is addressed
# txdb = makeTxDbFromGFF(annotation_gtf, format="gtf")
# txdb_exons = exonsBy(txdb, by="gene")
# txdb_exons = txdb_exons[names(txdb_exons) %in% rownames(dds)]
# rowRanges(dds) = sort(txdb_exons)
# export fpkm values
# fpkm = cbind.data.frame(, fpkm(dds, robust = F))
# fpkm = rownames_to_column(fpkm, 'gene')
# write.table(fpkm, file=snakemake@output[["fpkm"]], quote = F, row.names = F, sep = "\t")

# export deseq corrected sample expression values, merged with ids
sample_expression <- data.frame(DESeq2::counts(dds, normalized = TRUE),
  row.names = rownames(DESeq2::counts(dds))
)

sample_expression <- tibble::rownames_to_column(sample_expression, ids_in) %>%
  add_idcolnames(sp_arr)

write.table(
  sample_expression,
  file = snakemake@output[["sample_expression"]],
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

saveRDS(result_array, file = result_array_ids)
