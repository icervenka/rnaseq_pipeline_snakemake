#!/usr/bin/Rscript
suppressMessages(library(magrittr))
source("workflow/scripts/script_functions.R", local = TRUE)

rmarkdown::find_pandoc(dir = snakemake@params[["pandoc_path"]])

# snakemake parameters
experiment_name <- snakemake@params[["experiment_name"]]
fdr <- snakemake@params[["fdr"]]
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
species <- snakemake@params[["species"]]
ids_in <- snakemake@params[["gene_ids_in"]]
group <- snakemake@params[["mdplot_group"]]

dds <- readRDS(snakemake@input[[1]])
rld <- DESeq2::rlog(dds, blind = FALSE)
result_array <- readRDS(snakemake@input[[2]])
result_array_ids <- readRDS(snakemake@input[[3]])

report_layout <- paste0("_main_", snakemake@params[["report_layout"]])
outdir <- paste0("../../../", snakemake@params[["outdir"]])

knitr_output_options <- list(
  mathjax = NULL,
  self_contained = FALSE,
  lib_dir = paste0(outdir, "/libs")
)

# render html document
rmarkdown::render(
  paste0("workflow/scripts/deseq_report_rmd/", report_layout, ".Rmd"),
  output_file = paste0(snakemake@output[[1]]),
  output_dir = snakemake@params[["outdir"]],
  output_options = knitr_output_options,
  output_format = "all"
)

# # possibly move to separate file
# # Glimma MDS plot
# Glimma::glMDSPlot(dds,
#   path = snakemake@params[["outdir"]],
#   folder = "",
#   html = "mds-plot",
#   launch = FALSE
# )

# # possibly move to separate file
# # Glimma Volcano-Expression XY plots
# no_contrasts <- length(result_array)
# purrr::walk(1:no_contrasts, function(x) {
#   lfc_res <- result_array[[x]]
#   contrast <- names(result_array)[x]

#   lfc_res <- lfc_res[complete.cases(lfc_res), ]

#   counts <- counts(dds, normalized = TRUE)
#   counts <- counts[rownames(counts) %in% rownames(lfc_res), ]
#   counts <- counts[match(rownames(lfc_res), rownames(counts)), ]

#   annotation <- translate_gene_ids(rownames(lfc_res),
#     get_species_info(species),
#     from_type = type_in,
#     drop = FALSE
#   )
#   annotation <- annotation[match(rownames(lfc_res), annotation[[type_in]]), ]

#   status <- as.numeric(lfc_res$padj < fdr)
#   status <- status * sign(lfc_res$log2FoldChange)

#   Glimma::glXYPlot(
#     x = lfc_res$log2FoldChange,
#     # takes care of pvalues being 0 and transformed to Inf
#     y = -log10(lfc_res$pvalue + .Machine$double.xmin),
#     xlab = "log2FoldChange",
#     ylab = "-log10(pvalue)",
#     status = status,
#     counts = counts,
#     anno = data.frame(
#       GeneID = annotation$ENSEMBL,
#       Symbol = annotation$SYMBOL,
#       Name = annotation$GENENAME
#     ),
#     groups = dds[[group]],
#     samples = dds[["sample"]],
#     main = contrast,
#     side.main = "Symbol",
#     path = snakemake@params[["outdir"]],
#     folder = "",
#     html = paste0(contrast, "_expression"),
#     launch = FALSE
#   )
# })
