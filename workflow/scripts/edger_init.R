#!/usr/bin/Rscript
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(edgeR))

design_to_formula = function(design_string) {
  is_design_formula = startsWith(design_string, '~')
  if(is_design_formula) {
    f = as.formula(design_string)
  } else {
    f = as.formula(paste0("~ 0 + ", design_string))
  }
  return(f)
}

relevel_ref_levels = function(diffexp_object, levels) {
  ref_levels = str_split(levels, ',') %>% unlist
  if(length(ref_levels) > 1) {
    dobj = factor(diffexp_object, levels = ref_levels)
  } else {
    dobj = relevel(diffexp_object, ref_levels)
  }
  return(dobj)
}

# design = snakemake@params[["design"]]
# ref_levels = snakemake@params[["ref_levels"]]
# min_count = snakemake@params[["min_count"]]

setwd("~/rnaseq/lrt_hrt/")

design = "condition"
ref_levels = "LRT,LRTT,HRT,HRTT"
min_count = 10

# col_data = read.table(snakemake@input[["col_data"]], sep='\t', header = T)
col_data = read.table("metadata.tsv", sep='\t', header = T)
use_samples = paste(col_data$sample, collapse="|")
# count_data = read.table(snakemake@input[["count_data"]], sep='\t', header = T, row.names = 1)
count_data = read.table("counts/counts.tsv", sep='\t', header = T, row.names = 1)
count_data = count_data %>% dplyr::select(matches(use_samples))
names(count_data) = col_data$sample

targets = col_data %>% dplyr::select(-fq1, -fq2) %>% column_to_rownames(var = "sample")
mm = model.matrix(design_to_formula(design), data = targets)

# TODO if design is more complicated this might break
er = DGEList(counts=count_data)

er = er[filterByExpr(er), , keep.lib.sizes = F]
er = calcNormFactors(er)
er = estimateDisp(er, mm)

er$samples$group = relevel_ref_levels(er$samples$group, ref_levels)
# ref_levels = str_split(ref_levels, ',') %>% unlist
# if(length(ref_levels) > 1) {
#   er$samples$group <- factor(er$samples$group, levels = ref_levels)
# } else {
#   er$samples$group <- relevel(er$samples$group, ref_levels)
# }

group_levels = nlevels(col_data$condition)
level_comb = combn(levels(er$samples$group), 2)

design = model.matrix(~0+group, data=er$samples)
colnames(design) <- levels(er$samples$group)

if(group_levels > 1) {
  # glm approach
  glm_fit = glmQLFit(er, design)
  # save results of all contrasts in a list 
  result_array = apply(level_comb, 2, function(x) {
    contrast = makeContrasts(contrasts = paste0(x[2], "-", x[1]), levels=design)
    glmQLFTest(glm_fit, contrast=contrast)
  })
} else {
  # classic approach
  # save results of all contrasts in a list 
  result_array = apply(level_comb, 2, function(x) { 
    exactTest(er, pair=c(x[2], x[1]))
  })
}




