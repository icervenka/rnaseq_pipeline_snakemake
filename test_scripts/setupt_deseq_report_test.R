source("workflow/scripts/script_functions.R", local = TRUE)
library(magrittr)
library(yaml)

config = read_yaml("config.yaml")

experiment_name=config[["experiment_name"]]

fdr=config[['diffexp']][["fdr"]]
report_layout=config[["report"]][["layout"]]
individual=config[["report"]][["individual"]]
pca_groups=config[["report"]][["pca"]][["groups"]]
heatmap_mode=config[["report"]][["sample_heatmap"]][["mode"]]
heatmap_n_genes=config[["report"]][["sample_heatmap"]][["n_genes"]]
heatmap_cluster_metric_col=config[["report"]][["sample_heatmap"]][["cluster_metric_col"]]
heatmap_cluster_metric_row=config[["report"]][["sample_heatmap"]][["cluster_metric_row"]]
annotation_heatmap=config[["report"]][["sample_heatmap"]][["annotation"]]
annotation_sts=config[["report"]][["sample_distances_heatmap"]][["annotation"]]
upset_maxgroups=config[["report"]][["upset"]][["max_groups"]]
upset_min_group_size=config[["report"]][["upset"]][["min_group_size"]]
species=config[["species"]]
ids_in=config[['diffexp']][["input_gene_ids"]]
mdplot_group=config[["report"]][["mdplot_group"]]
pandoc_path=config[["pandoc_path"]]

outdir="diffexp/deseq_defaul/reports"

dds <- readRDS("diffexp/deseq_default/dds/dds.rds")
result_array <- readRDS("diffexp/deseq_default/dds/result_array.rds")
result_array_ids <- readRDS("diffexp/deseq_default/dds/result_array_ids.rds")

report_layout <- paste0("_main_", "flex")
outdir <- paste0("../../../", "diffexp/deseq_default/reports")

rmarkdown::find_pandoc(dir = "/home/ingwarr/miniconda3/envs/r_diffexp/bin")
# render html document
rmarkdown::render(
  paste0("workflow/scripts/deseq_report_rmd/", report_layout, ".Rmd"),
  output_file = "diffexp/deseq_default/reports/report.html",
  output_dir = "diffexp/deseq_default/reports/",
  output_options = knitr_output_options,
  output_format = "all"
)