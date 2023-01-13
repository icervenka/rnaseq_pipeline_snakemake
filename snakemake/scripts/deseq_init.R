#!/usr/bin/Rscript
suppressMessages(library(stringr))
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))

source("snakemake/scripts/script_functions.R", local = TRUE)

output <- snakemake@output[[1]]
# check if design is a formula, if not prepend a tilda
design <- modify_tilda(snakemake@params[["design"]])
ref_levels <- snakemake@params[["ref_levels"]]
min_count <- snakemake@params[["min_count"]]

col_data <- read.table(snakemake@input[["col_data"]], sep = "\t", header = TRUE)
col_data <- col_data %>%
  dplyr::select(-dplyr::any_of(fq, lane, read, sra)) %>%
  unique()
use_samples <- map_chr(
  col_data$sample,
  ~ paste0("\\.", .x, "\\.")
) %>% paste(collapse = "|")
count_data <- read.table(snakemake@input[["count_data"]],
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

dds <- dds[rowSums(counts(dds)) > min_count, ]

# keep for now
# remove_Y = snakemake@params[["remove_Y"]]
# if(remove_Y == TRUE) {
#   #todo for different genome assemblies
#   dds  = dds[ all(!seqnames(dds) %in% c("Y")), ]
# }

if (snakemake@params[["contrast_type"]] == "C" &&
  snakemake@params[["lfc_shrink"]] == "normal") {
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
  BPPARAM = MulticoreParam(snakemake@threads)
)

saveRDS(dds, file = output)
