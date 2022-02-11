# Bulk RNASeq processing (Snakemake pipeline)

This pipeline will process bulk RNASeq sequencing data from fastq files, bam
files or raw count matrices. Several most common aligning, read-counting programs
as well as programs for differential gene expression. In addition, tools for generating
QC reports, trimming and downloading SRA archives are implemented. Possibility of
defining simple custom workflows exists as well.

## Usage
You will need to copy or link fastq files in the `fastq` directory as well as provide
`metadata.tsv` with several required columns (See Input section).

Next the `config.yaml` file needs to be edited to reflect your desired output and
analysis parameters. Please see the parameter section of the readme file for detailed
description. `config_extra.yaml` file contains specific extra parameters that will be
passed to individual programs (aligners, read-counters etc.)

Run the snakemake pipeline by navigating to the folder with `Snakefile` and executing
it by command `snakemake`. One can optionally perform a test dry-run to see if pipeline
was correcly configured by running `snakemake -n`

Please see the documentation to [Snakemake pipeline language](https://snakemake.github.io/)
for detailed introduction how to run the program.

## Input
Fastq files with sequenced reads

`metadata.tsv` with several required columns:
- `sample` - sample name that represents one or more fastq files
- `fq` - name of fastq file, can be compressed
- `lane` - lane where the sample was one, multiple ones are allowed per sample
- `read` - '_R1' or '_R2' string at the end of the file basename in case the reads are paired
- `stranded` - indication whether stranded or unstranded protocols were used when
generating libraries. Values 'yes', 'no' and 'reverse' are accepted
- `group` - experimental group of a sample

Pipeline will merge files that are spread across multiple lanes and pair them correctly if necessary.

## Output
When run starting from fastq files, several output directories are created
- `align` contains aligned bam files
- `count` contains read counts
- `diffexp` contains data with differential expression divided into sub-analyses
- `logs` contains logs for all run processes as well as quality control reports

Note: Due to the fact that differential expression analysis usually needs more tailoring
based on experiment data, it is possible to add additional analyses without the need
to re-run alignment and read counting. It is sufficient to change the `diffexp: outdir`
parameter in the `config.yaml` file, change the parameters of DEG analysis such as
FDR or groups and it will be saved to its own `diffexp` subfolder. Parameters of the
analysis will be stored alongside.

## Configuration and Paramters

### Pipelines
Pipelines are defined in `snakemake/rules/pipelines.smk` file. New pipelines from
existing rules can be constructed if the output:input files of subsequent rules match.
Differential gene expression rules has to be defined as last. Skipping rules can be
achieved by specifying `skip_<process>.smk` in the pipelines.

Currently, following pipeline names are accepted, followed by applied rules
- cuffdiff: tophat, cufflinks, cuffdiff
- cuffdiff-denovo: star, cufflinks-denovo, cuffdiff
- stringtie: hisat, stringtie, ballgown
- deseq: star, featurecounts, deseq
- deseq-alt: hisat, featurecounts, deseq
- edger: star, featurecounts, edger
- limma: star, featurecounts, limma
- kallisto: kallisto, seluth
- only_download_sra: skips align, count and dge

### De-novo and guided assembly



## Other features

### Downloading SRA archives

### Creating coverage plots

### Trimming

### Quality control

## Issues

<!---
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
--->
