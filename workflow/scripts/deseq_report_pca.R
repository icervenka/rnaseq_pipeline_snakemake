#!/usr/bin/Rscript
suppressMessages(library(magrittr))
source("workflow/scripts/script_functions.R", local = TRUE)

# Parse snakemake items --------------------------------------------------------
## input
dds <- readRDS(snakemake@input[["dds"]])
## output
report_pca_out <- snakemake@output[["report_pca"]]
# params
outdir <- snakemake@params[["outdir"]]
dir_structure <- snakemake@params[["dir_structure"]]
sp_info <- get_species_info(snakemake@params[["species"]])
ids_in <- snakemake@params[["ids_in"]]
group <- snakemake@params[["group"]][1] # FIXME why [1]?
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

# FIXME
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

# set path for the main markdown file
report_layout_path <- paste0(
  "workflow/scripts/deseq_report_rmd/",
  "_pca_",
  "flex", # might change if other layouts are created later
  ".Rmd"
)

# output path is relative to this script that runs the markdown::render
## not snakemake base
relative_outdir <- paste0(dir_structure[["CD3UP"]], outdir)

knitr_output_options <- list(
  mathjax = NULL
)

rmarkdown::render(
  report_layout_path,
  output_file = paste0(dir_structure[["CD3UP"]], report_pca_out),
  output_options = knitr_output_options,
  output_format = "all"
)
