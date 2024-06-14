library(magrittr)
source("workflow/scripts/script_functions.R", local = TRUE)

# based on
# https://github.com/cole-trapnell-lab/cufflinks/issues/12
# for the count matrix to reflect raw counts cufflinks has to be most likely run
# with --no-length-correction, workflow taked care of it in case diffexp rules
# that require raw count matrices are specified, eg. DESeq, edgeR, limma.

# Parse snakemake items --------------------------------------------------------
## input
count_table <- snakemake@input[["count_table"]]
attr_table <- snakemake@input[["attr_table"]]
samples_table <- snakemake@input[["sample_table"]]
# count_table="counts/cuffnorm/genes.count_table"
# attr_table="counts/cuffnorm/genes.attr_table"
# samples_table="counts/cuffnorm/samples.table"
## output

## params
samples <- snakemake@params[["samples"]]
# samples = c("algwt", "lmcd1tg", "sra2", "algmut")

# Run --------------------------------------------------------------------------
# load the tables
counts <- read.delim(count_table)
attrs <- read.delim(attr_table)
cuffnorm_samples <- read.delim(samples_table)

# create unnormalization size factor vector
unnorm_vector <- cuffnorm_samples$internal_scale %>%
  setNames(cuffnorm_samples$sample_id)
unnorm_vector <- norm_vector[names(counts)[-1]]

# unnormalize the count matrix and round the counts to whole numbers
count_matrix <- round(sweep(counts[, -1], MARGIN = 2, unnorm_vector, `*`))

# reorganize and rename samples according to the provided parameter
sample_order <- purrr::map_int(samples, function(x) {
  stringr::str_which(names(count_matrix), x)
})

count_matrix <- count_matrix[sample_order]
names(count_matrix) <- samples

# create geneid column based on either gene symbols or cufflinks created genes
count_matrix_rownames <- attrs %>%
  dplyr::select(gene_id, gene_short_name) %>%
  dplyr::mutate(GENEID = ifelse(
    gene_short_name == "-",
    gene_id,
    gene_short_name
  ))

count_matrix <- cbind(GENEID = count_matrix_rownames$GENEID, count_matrix)

# Save -------------------------------------------------------------------------
write.table(
  count_matrix,
  snakemake@output[[1]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
