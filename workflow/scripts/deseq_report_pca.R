suppressMessages(library(magrittr))
source("workflow/scripts/script_functions.R", local = TRUE)

# Parse snakemake items --------------------------------------------------------
## input
dds <- readRDS(snakemake@input[["dds"]])
## output
report_pca_out = snakemake@input[["report_pca"]]
# params
sp_info = get_species_info(snakemake@params[["species"]])
ids_in <- snakemake@params[["ids_in"]]
group <- snakemake@params[["group"]][1] # TODO why [1]?
top_pcas <- snakemake@params[["top_pcas"]]
pca_ngenes <- snakemake@params[["ngenes"]]
top_gene_loadings <- snakemake@params[["top_gene_loadings"]]
pca2go_ngenes <- snakemake@params[["pca2go_ngenes"]]
pca2go_loadings_ngenes <- snakemake@params[["pca2go_loadings_ngenes"]]

# Run --------------------------------------------------------------------------
r <- DESeq2::rlog(dds, blind = TRUE)

pca_obj <- prcomp(t(SummarizedExperiment::assay(r)), rank. = 16)
pca_dims <- length(pca_obj$sdev)
top_pcas <- min(top_pcas, pca_dims)

pca_corr <- pcaExplorer::correlatePCs(
  pca_obj,
  SummarizedExperiment::colData(dds),
  pcs = 1:pca_dims
)

# TODO fixme
# anno_df_orgdb <- get_annotation_orgdb(
#   dds = dds,
#   orgdb_species = sp_info["orgdb"],
#   idtype = id_type
# )

pca2go <- pcaExplorer::limmaquickpca2go(
  r,
  pca_ngenes = pca2go_ngenes,
  loadings_ngenes = pca2go_loadings_ngenes,
  inputType = ids_in,
  organism = sp_info["orgdb_short"]
)


knitr_output_options <- list(
  mathjax = NULL
)

rmarkdown::render(paste0("workflow/scripts/deseq_report_rmd/_pca_flex.Rmd"),
  output_file = paste0(CD3UP, report_pca_out),
  output_options = knitr_output_options,
  output_format = "all"
)
