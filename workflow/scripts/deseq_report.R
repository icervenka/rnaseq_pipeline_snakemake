#!/usr/bin/Rscript
suppressMessages(library(magrittr))
source("workflow/scripts/script_functions.R", local = TRUE)

# Parse snakemake items --------------------------------------------------------
## input
dds <- readRDS(snakemake@input[["dds"]])
result_array <- readRDS(snakemake@input[["result_array"]])
result_array_ids <- readRDS(snakemake@input[["result_array_ids"]])
## output
report_out <- snakemake@output[["report"]]
## params
experiment_name <- snakemake@params[["experiment_name"]]
outdir <- snakemake@params[["outdir"]]
species <- snakemake@params[["species"]]
ids_in <- snakemake@params[["ids_in"]]
fdr <- snakemake@params[["fdr"]]
report_layout <- snakemake@params[["report_layout"]]
individual <- snakemake@params[["individual"]]
pca_groups <- snakemake@params[["pca_groups"]]
heatmap_mode <- snakemake@params[["heatmap_mode"]]
heatmap_n_genes <- snakemake@params[["heatmap_n_genes"]]
heatmap_cluster_metric_col <- snakemake@params[["heatmap_cluster_metric_col"]]
heatmap_cluster_metric_row <- snakemake@params[["heatmap_cluster_metric_row"]]
annotation_heatmap <- snakemake@params[["annotation_heatmap"]]
annotation_sts <- snakemake@params[["annotation_sts"]]
upset_maxgroups <- snakemake@params[["upset_maxgroups"]]
upset_min_group_size <- snakemake@params[["upset_min_group_size"]]

# Run --------------------------------------------------------------------------
rld <- DESeq2::rlog(dds, blind = FALSE)

# set path for the main markdown file
report_layout_path <- paste0(
  "workflow/scripts/deseq_report_rmd/",
  "_main_",
  report_layout,
  ".Rmd"
)

# output path is relative to this script that runs the markdown::render
## not snakemake base
relative_outdir <- paste0(CD3UP, outdir)

knitr_output_options <- list(
  mathjax = NULL,
  self_contained = FALSE,
  lib_dir = paste0(relative_outdir, "/libs")
)

# FIXME gives warning for diffexp dir creation permission denied, but works ok
# render html document
rmarkdown::render(
  report_layout_path,
  output_file = paste0(CD3UP, report_out),
  output_options = knitr_output_options,
  output_format = "all"
)
