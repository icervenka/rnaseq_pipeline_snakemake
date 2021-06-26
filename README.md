# align options ----------------------------------------------------------------
# sra archive download
# requires that the column sra with SRA accession is present in metadata file
# names of the files in fq column have to follow the scheme:
# {name}{read}.fastq.gz, name is the SRA accession and
# read is empty string for unpaired data and _1 and _2 for paired mates
# some filesystems don't support zcat pipe into eg STAR, alignment will fail

## pipelines are defined in snakemake/rules/pipelines.smk
## following pipeline names are accepted, followed by applied rules
# cuffdiff : tophat, cufflinks, cuffdiff
# cuffdiff-denovo: star, cufflinks-denovo, cuffdiff
# stringtie : hisat, stringtie, ballgown
# deseq : star, featurecounts, deseq
# deseq-alt : hisat, featurecounts, deseq
# edger : star, featurecounts, edger
# limma : star, featurecounts, limma
# kallisto : kallisto, seluth
# only_download_sra : skips align, count and deg

# some specific config parameters such as extra parameters to aligners
# or diffexp analysis are present in snakemake/scripts/extra_config.yaml

# specify which comparisons to include, can be one of:
# - type 'A' - specify list of columns in metadata file, all combinations of
#   contrasts will be generated
#   eg. ["condition", "gender"]
# - type 'B' - named lists of three items, name of grouping and two ref_level
#   combinations. The baseline is specified as second. Name has to be unique
#   eg. treatment: ["condition", "treat", "ctrl"]
# - type 'C' - lists of contrasts to include, has to correspond
#   to resultsNames of dds object. Names have to be unique
#   eg. comparison1: ['condition_Trt_vs_Ctrl', 'genotypeMU.conditionTrt']
# - type 'C' - numeric lists of contrasts to include and their weights,
#   has to correspond to resultsNames of dds object. Names have to be unique
#   eg. comparison1: [-1, 1, 0, 0, 1]
# - type 'D' - full and reduced model formula in a list (sleuth only)
# - type 'E' - link to file with design matrix
