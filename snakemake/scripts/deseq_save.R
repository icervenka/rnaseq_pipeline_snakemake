#!/usr/bin/Rscript

suppressMessages(library(tibble))
suppressMessages(library(DESeq2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(dplyr))

source("snakemake/scripts/script_functions.R", local = TRUE)

sp_arr <- get_species_info(snakemake@params[["species"]])
outdir <- snakemake@params[["outdir"]]
ids_in <- snakemake@params[["gene_ids_in"]]
result_array <- readRDS(snakemake@input[["result_array"]])
dds <- readRDS(snakemake@input[["dds"]])

# TODO validate translate id function
# create a dataframe with IDs
ids <- suppressMessages(
  clusterProfiler::bitr(
    rownames(dds),
    fromType = ids_in,
    toType = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"),
    OrgDb = sp_arr["orgdb"]
  )
)

result_array <- lapply(result_array, function(x) {
  x <- tibble::rownames_to_column(data.frame(x), ids_in)
  x <- merge(ids, x, by = ids_in)
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
    result_array[[x]] %>% filter(padj < 0.05),
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
sample_expression <- data.frame(counts(dds, normalized = TRUE),
  row.names = rownames(counts(dds))
)
sample_expression <- rownames_to_column(sample_expression, "gene")
sample_expression <- merge(ids, sample_expression, by.x = ids_in, by.y = "gene")
write.table(
  sample_expression,
  file = snakemake@output[["sample_expression"]],
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

saveRDS(result_array, file = snakemake@output[["result_array_ids"]])
