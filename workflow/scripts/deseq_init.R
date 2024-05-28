#!/usr/bin/Rscript
suppressMessages(library(magrittr))
source("workflow/scripts/script_functions.R", local = TRUE)

# Parse snakemake items --------------------------------------------------------
## inputs
col_data <- snakemake@input[["col_data"]]
count_data <- snakemake@input[["count_data"]]
## outputs
output <- snakemake@output[["dds"]]
## params
# check if design is a formula, if not prepend a tilda
design <- modify_tilda(snakemake@params[["design"]])
ref_levels <- snakemake@params[["ref_levels"]]
min_count <- snakemake@params[["min_count"]]
contrast_type <- snakemake@params[["contrast_type"]]
lfc_shrink <- snakemake@params[["lfc_shrink"]]
threads <- snakemake@threads

# Run --------------------------------------------------------------------------
col_data <- read.table(col_data, sep = "\t", header = TRUE)
col_data <- col_data %>%
  dplyr::select(-dplyr::any_of(c("fq", "lane", "read", "sra"))) %>%
  unique()

use_samples <- purrr::map_chr(
  col_data$sample,
  ~ paste0("\\.", .x, "\\.")
) %>%
  paste(collapse = "|")

count_data <- read.table(
  count_data,
  sep = "\t",
  header = TRUE,
  row.names = 1
)

count_data <- count_data %>% dplyr::select(matches(use_samples))
names(count_data) <- unique(col_data$sample)

dds <- DESeq2::DESeqDataSetFromMatrix(
  count_data,
  col_data,
  design = as.formula(design)
)

# relevel the factors for dds object
# set reference level if only one is supplied, else relevel the whole factor
for (x in names(ref_levels)) {
  if (length(ref_levels[[x]]) > 1) {
    dds[[x]] <- factor(dds[[x]], levels = ref_levels[[x]])
  } else {
    dds[[x]] <- relevel(dds[[x]], ref_levels[[x]])
  }
}

dds <- dds[rowSums(DESeq2::counts(dds)) > min_count, ]

# keep for now
# remove_Y = snakemake@params[["remove_Y"]]
# if(remove_Y == TRUE) {
#   #todo for different genome assemblies
#   dds  = dds[ all(!seqnames(dds) %in% c("Y")), ]
# }

if (contrast_type == "C" && lfc_shrink == "normal") {
  set_beta_prior <- TRUE
  warning("You have selected numeric contrasts with 'normal'
  shrinkage. This runs DESeq with betaPrior = T.")
} else {
  set_beta_prior <- FALSE
}

dds <- DESeq2::DESeq(dds,
  betaPrior = set_beta_prior,
  quiet = TRUE,
  parallel = TRUE,
  BPPARAM = BiocParallel::MulticoreParam(threads)
)

# Save -------------------------------------------------------------------------
saveRDS(dds, file = output)
