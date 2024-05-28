source("workflow/scripts/script_functions.R", local = TRUE)
library(magrittr)

result_array <- readRDS("diffexp/deseq_default/dds/result_array.rds")
dds <- readRDS("diffexp/deseq_default/dds/dds.rds")
outdir <- "diffexp/deseq_default"
result_array_ids <- "diffexp/deseq_default/dds/result_array_ids.rds"
sample_expression <- "diffexp/deseq_default/sample_expression.csv"
sp_info <- "mouse"
ids_in <- "ENSEMBL"
