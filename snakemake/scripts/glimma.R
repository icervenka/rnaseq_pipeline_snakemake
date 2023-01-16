library(Glimma)
library(clusterProfiler)
source("snakemake/scripts/script_functions.R")





dds = dds
result_array = result_array
no_contrasts = length(result_array)

glMDSPlot(dds,
          path = "",
          folder = "",
          html = "")


species = "mouse" # to param
group = "group"
fdr = 0.05
type_in = "ENSEMBL"

x = 1
walk(1:no_contrasts, function(x) {
  lfc_res = result_array[[x]]
  contrast = names(result_array)[x]

  lfc_res = lfc_res[complete.cases(lfc_res),]
  
  counts = counts(dds, normalized = TRUE)
  counts = counts[rownames(counts) %in% rownames(lfc_res), ]
  counts = counts[match(rownames(lfc_res), rownames(counts)), ]
  
  annotation = translate_gene_ids(rownames(lfc_res),
                                  get_species_info(species),
                                  from_type = type_in,
                                  drop = F)
  annotation = annotation[match(rownames(lfc_res), annotation[[type_in]]),]
  
  status = as.numeric(lfc_res$padj < fdr)
  status = status * sign(lfc_res$log2FoldChange)
  
  glXYPlot(x = lfc_res$log2FoldChange, 
           y = -log10(lfc_res$pvalue),
           xlab = "log2FoldChange",
           ylab = "-log10(pvalue)",
           status = status,
           counts = counts,
           anno = data.frame(GeneID = annotation$ENSEMBL, 
                             Symbol = annotation$SYMBOL, 
                             Name = annotation$GENENAME), 
           side.main = "Symbol",
           groups = dds$group,
           samples = dds$sample,
           main = contrast
           # path = "",
           # folder = "",
           # html = ""
  )
  
})