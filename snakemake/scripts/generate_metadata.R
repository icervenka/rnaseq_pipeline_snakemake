#!/usr/bin/Rscript
library(stringr)
library(tibble)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 6) {
  stop("Exactly 6 arguments must be supplied")
}

in.path = args[1]
out.path = paste0(args[2], "/")
pattern = args[3]
remove_extension = args[4]
into = trimws(unlist(strsplit(args[5], ",")))
into[into == "NA"] = NA
sep = args[6]

files = list.files(path = in.path, pattern = pattern) %>%
  tibble(fq = .)

basenames = str_replace_all(files$fq, remove_extension, "") %>%
  tibble(fq = .) %>%
  separate(fq, into = into, sep = sep, remove = T)

metadata = cbind(files, basenames)
write.table(metadata, paste0(out.path, "metadata.tsv"), sep = "\t", quote = F, row.names = F)
