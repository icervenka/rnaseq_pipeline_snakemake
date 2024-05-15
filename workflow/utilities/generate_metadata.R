#!/usr/bin/Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))

option_list = list(
  make_option(c("-o", "--output"), type="character", default="metadata.tsv", 
              help="output file to which metadata will be written [default=%default]", 
              metavar="character"),
  make_option(c("--ext"), type="character", default=".fastq.gz", 
              help="extension of the files that will be removed to generate base
              names [default=%default]", 
              metavar="character"),
  make_option(c("--columns"), type="character", default="NA,sample,NA,NA,lane,read", 
              help="comma separated name of descriptor columns created from a filename. 
              Use NA to skip the column. Number of descriptors has to be equal to 
              number of fields created by separating at --sep. Sample descriptor
              is required [default=%default]", 
              metavar="character"),
  make_option(c("--sep"), type="character", default="_", 
              help="character separator used to split filename into columns [default=%default]", 
              metavar="character")
); 

opt_parser = OptionParser(usage = "%prog [options] input_directory", option_list=option_list);
opt = parse_args(opt_parser, positional_arguments = 1);

opt$options$columns = trimws(unlist(strsplit(opt$options$columns, ",")))
opt$options$columns[opt$options$columns == "NA"] = NA

files = list.files(path = opt$args) %>%
  tibble::tibble(fq = .)

metadata = files %>% 
  dplyr::mutate(fq_basename = stringr::str_replace_all(fq, opt$options$ext, "")) %>%
  tidyr::separate(fq_basename, into = opt$options$columns, sep = opt$options$sep, remove = T) %>%
  dplyr::select(sample, fq, everything())

if(is.null(metadata[['group']])) {
  metadata = metadata %>%
    dplyr::mutate(group = sample)
}

write.table(metadata, 
            opt$options$o, 
            sep = "\t", 
            quote = F, 
            row.names = F)
