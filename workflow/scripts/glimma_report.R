#!/usr/bin/Rscript
library(magrittr)
source("workflow/scripts/script_functions.R", local = TRUE)

# Parse snakemake items --------------------------------------------------------
## input
dds <- readRDS(snakemake@input[["dds"]])
result_array <- readRDS(snakemake@input[["result_array"]])
no_contrasts <- length(result_array)
## output

## params
outdir <- snakemake@params[["outdir"]]
species <- snakemake@params[["species"]]
ids_in <- snakemake@params[["diffexp"]][["ids_in"]]
fdr <- snakemake@params[["diffexp"]][["fdr"]]
group <- snakemake@params[["report"]][["mdplot_group"]]

diffexp_extra <- snakemake@params[["diffexp_extra"]]
report_extra <- snakemake@params[["report_extra"]]
# Run --------------------------------------------------------------------------
Glimma::glMDSPlot(
  dds,
  path = outdir,
  folder = "",
  html = "mds-plot",
  launch = FALSE
)

purrr::walk(1:no_contrasts, function(x) {
  lfc_res <- result_array[[x]] %>%
    as.data.frame()
  contrast <- names(result_array)[x]

  lfc_res <- lfc_res[complete.cases(lfc_res), ]

  cnts <- DESeq2::counts(dds, normalized = TRUE)
  cnts <- cnts[rownames(cnts) %in% rownames(lfc_res), ]
  cnts <- cnts[match(rownames(lfc_res), rownames(cnts)), ]

  annotation <- translate_gene_ids(rownames(lfc_res),
    get_species_info(species),
    from_type = ids_in,
    drop = FALSE
  )

  annotation <- annotation[match(rownames(lfc_res), annotation[[ids_in]]), ]

  status <- as.numeric(lfc_res$padj < fdr)
  status <- status * sign(lfc_res$log2FoldChange)

  Glimma::glXYPlot(
    x = lfc_res$log2FoldChange,
    y = -log10(lfc_res$pvalue),
    xlab = "log2FoldChange",
    ylab = "-log10(pvalue)",
    status = status,
    counts = cnts,
    anno = data.frame(
      GeneID = annotation$ENSEMBL,
      Symbol = annotation$SYMBOL,
      Name = annotation$GENENAME
    ),
    side.main = "Symbol",
    groups = dds$group,
    samples = dds$sample,
    main = contrast,
    path = outdir,
    folder = "",
    html = paste0(contrast, "_expression"),
    launch = FALSE
  )
})
