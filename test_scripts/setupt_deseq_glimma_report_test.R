source("workflow/scripts/script_functions.R", local = TRUE)
library(magrittr)
library(yaml)

config <- read_yaml("config.yaml")

species=config[['species']]
id_type=config[['diffexp']][["input_gene_ids"]]
group=config[['report']][["pca"]][['groups']]
ngenes=config[['report']][["pca"]][['ngenes']]
top_pcas=config[['report']][["pca"]][['top_pcas']]
top_gene_loadings=config[['report']][["pca"]][['top_gene_loadings']]
pca2go_ngenes=config[['report']][["pca"]][['pca2go_ngenes']]
pca2go_loadings_ngenes=config[['report']][["pca"]][['pca2go_loadings_ngenes']]
pandoc_path=config[["pandoc_path"]]

dds <- dds <- readRDS("diffexp/deseq_default/dds/dds.rds")
sp_info <- get_species_info(species)

r <- DESeq2::rlog(dds, blind = TRUE)
pca_obj <- prcomp(t(SummarizedExperiment::assay(r)), rank. = 16)
pca_dims <- length(pca_obj$sdev)
pca_corr <- pcaExplorer::correlatePCs(pca_obj, SummarizedExperiment::colData(dds), pcs = 1:pca_dims)


top_pcas <- min(top_pcas, pca_dims)
plts <- purrr::map(combn(1:top_pcas, 2) %>% as.data.frame(), function(x) {
  pcaExplorer::pcaplot(
    r,
    intgroup = group,
    ntop = ngenes,
    text_labels = FALSE,
    pcX = x[1],
    pcY = x[2]
  ) +
    ggplot2::theme(legend.position = "none")
})

  pcaExplorer::pcaplot(
    r,
    intgroup = group,
    ntop = ngenes,
    text_labels = FALSE,
    pcX = 4,
    pcY = 5
  )