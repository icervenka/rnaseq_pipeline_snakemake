source("workflow/scripts/script_functions.R", local = TRUE)

col_data <- "metadata.tsv"
count_data <- "counts/count_matrix.txt"
output <- ""
# check if design is a formula, if not prepend a tilda
design <- modify_tilda("group")
ref_levels <- c("alg, lmcd, liver")
min_count <- 10
contrast_type <- "A"
lfc_shrink <- "normal"
threads <- 1

snakemake <- readRDS(file = "snakemake.rds")



unpack_named_snakemake_slots <- function(snakemake, slots = c("input", "output", "params")) {
    lapply(slots, function(x) {
        named_vals = vapply(names(attributes(snakemake)[[x]]), nchar, FUN.VALUE = numeric(1)) > 0
        lapply(names(attributes(snakemake)[[x]][named_vals]), function(y) {
            assign(y, attributes(snakemake)[[x]][[y]])
        })
    })
}


unpack_snakemake_slots(snakemake)
