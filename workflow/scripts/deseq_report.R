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
dir_structure <- snakemake@params[["dir_structure"]]
outdir <- snakemake@params[["outdir"]]
species <- snakemake@params[["species"]]

ids_in <- snakemake@params[["diffexp"]][["ids_in"]]
fdr <- snakemake@params[["diffexp"]][["fdr"]]

report_params <- snakemake@params[["report"]]
report_layout <- report_params[["report_layout"]]
individual <- report_params[["individual"]]
self_contained <- report_params[["self_contained"]]
pca_groups <- report_params[["pca_groups"]]
heatmap_mode <- report_params[["sample_heatmap"]][["mode"]]
heatmap_n_genes <- report_params[["sample_heatmap"]][["ngenes"]]
heatmap_cluster_metric_col <- report_params[["sample_heatmap"]][["cluster_metric_col"]]
heatmap_cluster_metric_row <- report_params[["sample_heatmap"]][["cluster_metric_row"]]
annotation_heatmap <- report_params[["sample_heatmap"]][["annotation"]]
annotation_sts <- report_params[["sample_distance_heatmap"]][["annotation"]]
upset_maxgroups <- report_params[["upset"]][["maxgroups"]]
upset_min_group_size <- report_params[["upset"]][["min_group_size"]]

diffexp_extra <- snakemake@params[["diffexp_extra"]]
report_extra <- snakemake@params[["report_extra"]]
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
relative_outdir <- paste0(dir_structure[["CD3UP"]], outdir)

knitr_output_options <- list(
  mathjax = NULL,
  self_contained = self_contained, # FALSE,
  lib_dir = paste0(relative_outdir, "/libs")
)

rmarkdown::render(
  report_layout_path,
  output_file = paste0(dir_structure[["CD3UP"]], report_out),
  output_options = knitr_output_options,
  output_format = "all"
)
