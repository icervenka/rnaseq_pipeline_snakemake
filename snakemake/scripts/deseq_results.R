#!/usr/bin/Rscript
suppressMessages(library(DESeq2))
suppressMessages(library(tibble))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(tidyverse))

source("snakemake/scripts/common.R")

# TODO add parallelization to lfcshrink
dds = readRDS(snakemake@input[["dds"]])
lfc_shrink = snakemake@params[["lfc_shrink"]]
fdr = snakemake@params[["fdr"]]
contrasts = snakemake@params[["contrasts"]]
contrast_names = names(contrasts)
contrast_type = snakemake@params[["contrast_type"]]

if(contrast_type == 'A') {
  sc = contrasts %>% unlist
  selected_contrasts = map_dfr(sc, function(x) {
    combn(levels(dds[[x]]), 2) %>% 
      t %>% 
      as.data.frame(stringsAsFactors = F) %>%
      dplyr::mutate(cond = x) %>%
      dplyr::select(cond, V2, V1)
  }) 
  contrast_names = selected_contrasts %>%
    dplyr::mutate(name = name_contrast(V2, V1)) %>%
    pull(name)
  selected_contrasts = selected_contrasts %>%
    t %>% 
    as.data.frame(stringsAsFactors = F)
} else if(contrast_type == 'B' | contrast_type == 'C') {
  selected_contrasts = dplyr::bind_cols(contrasts)
} else {
  stop("Incorrect contrast type specified.")
}

print(selected_contrasts)

result_array = map(selected_contrasts, function(x) { 
  sc = x 
  res = results(dds, contrast=sc, alpha = fdr)
  if(lfc_shrink == "apeglm") {
    suppressMessages(library(apeglm))
    #TODO
    lfc_res = lfcShrink(dds, contrast=sc, res=res, type = "ashr")
  } else if(lfc_shrink == "ashr") {
    suppressMessages(library(ashr))
    lfc_res = lfcShrink(dds, contrast=sc, res=res, type = "ashr")
  } else if(lfc_shrink == "normal" & contrast_type == "C"){
    message("When selecting numeric contrast and normal lfc shrinkage need following number of contrasts specified:")
    print(resultsNames(dds))
    lfc_res = res
  } else {
    lfc_res = lfcShrink(dds, contrast=sc, res=res, type = "normal")
  }
  lfc_res[order(lfc_res$padj),]
})

names(result_array) = contrast_names
saveRDS(result_array, file=snakemake@output[[1]])

# save contrasts to a file
# write(names(contrasts), snakemake@output[[2]])
