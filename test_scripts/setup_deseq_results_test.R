source("workflow/scripts/script_functions.R", local = TRUE)

dds <- readRDS("diffexp/deseq_default/dds/dds.rds")
lfc_shrink <- "normal"
fdr <- 0.05
contrasts <- "group"
contrast_names <- names(contrasts)
contrast_type <- "A"
threads <- 1

