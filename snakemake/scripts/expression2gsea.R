dds_to_gsea_norm = function(dds, group, filename) {
  cnts = counts(dds, normalized = T) %>%
    as_tibble(rownames = "NAME") %>%
    dplyr::mutate(DESCRIPTION = "NA") %>%
    dplyr::select(NAME, DESCRIPTION, everything())

  no_classes = colData(dds)[[group]] %>% unique %>% length
  class_labels = colData(dds)[[group]] %>% unique %>% as.character %>% rev
  no_samples = length(colData(dds)$sample)
  sample_labels = colData(dds)[[group]] %>% as.character %>% rev

  sink(paste0(filename, ".cls"))
  cat(paste(no_samples, no_classes, 1, sep = " "))
  cat('\n')
  cat(paste("#", paste(class_labels, collapse = " "), sep = " "))
  cat('\n')
  cat(paste(sample_labels, collapse = " "))
  cat('\n')
  sink()

  write.table(cnts, paste0(filename, ".txt"), quote = F, sep = '\t', row.names = F)
}
