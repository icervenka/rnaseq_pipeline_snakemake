source("workflow/scripts/script_functions.R", local = TRUE)
library(magrittr)
library(yaml)

config <- read_yaml("config.yaml")

outdir="diffexp/deseq_default/reports/"
species=config[['species']]
ids_in=config[['diffexp']][["input_gene_ids"]]
group=config[["report"]][["mdplot_group"]]
fdr=config[['diffexp']][["fdr"]]

dds <- readRDS("diffexp/deseq_default/dds/dds.rds")
result_array <- readRDS("diffexp/deseq_default/dds/result_array.rds")
no_contrasts <- length(result_array)
sp_info <- get_species_info(species)

Glimma::glMDSPlot(
  dds,
  path = outdir,
  folder = "",
  html = "mds-plot",
  launch = FALSE
)

purrr::walk(1:no_contrasts, function(x) {
  lfc_res <- result_array[[x]]
  contrast <- names(result_array)[x]

  lfc_res <- lfc_res[complete.cases(lfc_res), ]

  counts <- counts(dds, normalized = TRUE)
  counts <- counts[rownames(counts) %in% rownames(lfc_res), ]
  counts <- counts[match(rownames(lfc_res), rownames(counts)), ]

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
    counts = counts,
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
