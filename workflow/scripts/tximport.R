suppressMessages(library(magrittr))
source("workflow/scripts/script_functions.R", local = TRUE)

# Parse snakemake items --------------------------------------------------------
## input
files <- snakemake@input[["files"]]
gtf <- snakemake@input[["gtf"]]
## output

## params
metadata <- read.table(
  snakemake@params[["metadata"]],
  sep = "\t",
  header = TRUE
)
samples <- snakemake@params[["samples"]]
sp_arr <- get_species_info(snakemake@params[["species"]])
pipeline <- snakemake@params[["pipeline"]]
tool <- snakemake@params[["tool"]]
ids_in <- snakemake@params[["ids_in"]]
out_type <- snakemake@params[["out_type"]]
# TODO implement extra parameters
# extra =

# Run --------------------------------------------------------------------------
files <- setNames(files, samples)

txdb <- GenomicFeatures::makeTxDbFromGFF(
  file = gtf,
  format = "gtf",
  taxonomyId = as.numeric(sp_arr[["tax"]])
)

k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

tx2gene <- tx2gene %>%
  dplyr::mutate(TXNAME = remove_ens_gene_version(TXNAME)) %>%
  dplyr::mutate(GENEID = remove_ens_gene_version(GENEID))

if (tool == "" || tool == "auto") {
  tximport_type <- set_tximport_type(pipeline)
} else {
  tximport_type <- tool
}

# if one has quantified with Salmon or alevin with human or mouse transcriptomes
# is to use the tximeta function from the tximeta Bioconductor package
txi <- tximport::tximport(files,
  type = tximport_type,
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE,
  ignoreAfterBar = TRUE
  # !!!extra
)

count_matrix <- txi[[out_type]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = ids_in)

# Save -------------------------------------------------------------------------
write.table(
  count_matrix,
  snakemake@output[[1]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
