suppressMessages(library(rmarkdown))
suppressMessages(library(pcaExplorer))
suppressMessages(library(DESeq2))
suppressMessages(library(plotly))
suppressMessages(library(heatmaply))
suppressMessages(library(DT))
suppressMessages(library(tidyverse))
source("snakemake/scripts/common.R", local = TRUE)

# snakemake params ----------------------------------------------------------
species <- snakemake@params[["species"]]
id_type <- snakemake@params[["id_type"]]
group <- snakemake@params[["group"]][1]
top_pcas <- snakemake@params[["top_pcas"]]
pca_ngenes <- snakemake@params[["ngenes"]]
top_gene_loadings <- snakemake@params[["top_gene_loadings"]]
pca2go_ngenes <- snakemake@params[["pca2go_ngenes"]]
pca2go_loadings_ngenes <- snakemake@params[["pca2go_loadings_ngenes"]]

# snakemake inputs ----------------------------------------------------------
dds <- readRDS(snakemake@input[["dds"]])
sp_info <- get_species_info(species)

# prepare files -------------------------------------------------------------
r <- DESeq2::rlog(dds, blind = TRUE)
pca_obj <- prcomp(t(assay(r)), rank. = 16)
pca_dims <- length(pca_obj$sdev)
pca_corr <- pcaExplorer::correlatePCs(pca_obj, colData(dds), pcs = 1:pca_dims)
anno_df_orgdb <- get_annotation_orgdb(
  dds = dds,
  orgdb_species = sp_info["orgdb"],
  idtype = id_type
)

pca2go <- pcaExplorer::limmaquickpca2go(
  r,
  pca_ngenes = pca2go_ngenes,
  loadings_ngenes = pca2go_loadings_ngenes,
  inputType = id_type,
  organism = sp_info["orgdb_short"]
)

# run pca report ------------------------------------------------------------
knitr_output_options <- list(
  mathjax = NULL
)

render(paste0("snakemake/scripts/deseq_report_rmd/_pca_flex.Rmd"),
  output_file = paste0("../../../", snakemake@output[[1]]),
  output_options = knitr_output_options,
  output_format = "all"
)
