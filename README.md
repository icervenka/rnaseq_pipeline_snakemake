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

Currently, following pipeline names are accepted, followed by applied rules.
 - `download_only` - skip_align, skip_count, skip_diffexp
 - `star_only` - star, skip_count, skip_diffexp
 - `hisat_only` - hisat, skip_count, skip_diffexp
 - `tophat_only` - tophat, skip_count, skip_diffexp
 - `bowtie2_only` - bowtie2, skip_count, skip_diffexp
 - `salmon_only` - salmon, skip_count, skip_diffexp
 - `kallisto_only` - kallisto, skip_count, skip_diffexp
 - `featurecounts_only` - skip_align, featurecounts, skip_diffexp
 - `htseq_only` - skip_align, htseq, skip_diffexp
 - `stringtie_only` - skip_align, stringtie, skip_diffexp
 - `stringtie_expression_only` - skip_align, stringtie_expression, skip_diffexp
 - `cufflinks_only` - skip_align, cufflinks, skip_diffexp
 - `deseq_only` - skip_align, featurecounts, deseq
 - `deseq` - star, featurecounts, deseq

### De-novo and guided assembly

### config.yaml parameters

#### general
- `experiment_name` - name of the experiment (no whitespaces)
- `metadata` - names of the file where experiment metadata is specified
- `threads` - maximum amount of threads the pipeline is allowed to use
- `species` - source species of the experiment. This is mainly used during conversion of gene IDs, human, mouse and rat are currently supported.
- `pipeline` - name of the pipeline to use for the analysis. Can be any that is specified in the `snakemake/rules/pipelines.smk` file
- `index` - location of genome index for selected aligner
- `gtf` - location of gtf file with features to be counted, chromosome naming and locations have to be compatible with the index used
- `fasta` - fasta file used to generate the index. This is required for some aligners such as bowtie2

#### trim

- `adapters_single` - see [cutadapt adapter types](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types) for details
- `adapters_paired` - see [cutadapt adapter types](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types) for details
- `quality` - see [cutadapt quality trimming](https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming) for details
- `extra` - extra parameters to cutadapt


####  coverage
- `calulate` - whether coverage should be calculated. accepts "yes" or "no" string. This solution is a bit cumbersome, but works for now
- `split_strands` -
- `bin_size` -
- `normalize_using` - methods for normalizing reads,
- `extra` - extra parameters to coverage calculations

#### count
- `multimap` - if multimapping reads should be counted
- `overlap` - if reads belonging to the overlapping genes should be counted. If more granular control is needed about how to count overlaps, these need to be specified in `config_extra.yaml`

#### diffexp

- `outdir` - name of the analysis and final output directory. If multiple analyses of differential expression need to be conducted without repeating the alignment and read counting, it is sufficient to change the name of this directory along with the analysis parameters and it will be added to `diffexp` directory
- `design` - specify the experimental design. the default value for one design variable is 'group'. If only one grouping variable is present, specify just the name of the column in `metadata.tsv` file. If you wish to use formula to combine multiple groups, mark the beginning with ~ any use combination of columns present in `metadata.tsv` file.
- `reference_levels` - ordering of the reference levels in the `design` column from which contrasts will be constructed. If only one level is specified it is taken as baseline. If empty, alphabetical order will be used.
- `contrast_type` - specify which type of comparisons to include, can be one of:
  - type 'A' - specify list of columns in metadata file, all combinations of
contrasts will be generated
eg. ["condition", "gender"]
  - type 'B' - named lists of three items, name of grouping and two ref_level
combinations. The baseline is specified as second. Name has to be unique
eg. treatment: ["condition", "treat", "ctrl"]
  - type 'C' - lists of contrasts to include, has to correspond
to resultsNames of dds object. Names have to be unique (DESeq2 only).
eg. comparison1: ['condition_Trt_vs_Ctrl', 'genotypeMU.conditionTrt']
  - type 'C' - numeric lists of contrasts to include and their weights,
has to correspond to resultsNames of dds object. Names have to be unique
eg. comparison1: [-1, 1, 0, 0, 1]
  - type 'D' - full and reduced model formula in a list (sleuth only)
  - type 'E' - link to file with design matrix
- `contrasts` - list of actual contrast, using naming conventions specified above
- `fdr` - FDR cutoff to consider for diferentially expressed genes
- `gene_min_readcount` - filter out genes that have smaller total amount of reads across all samples
- `input_gene_ids` - type of gene IDs present in the input file (generally the read counts). Can be one of supported by org.db, usually "ENSEMBL", "ENTREZ", "SYMBOL".
- `output_gene_ids` - type of gene IDs to use in output files.
- `drop_no_sybmol` - drop genes from output that don't map to any gene symbol, otherwise they will be replaced by gene ID.




### config_extra.yaml parameters
Contains some specific config parameters such as extra parameters to aligners, read counting etc.



## Other features

### Downloading SRA archives
Rule that will download file from SRA archive is available. In order to download SRR files, a `metadata.tsv` file has to constructed as follows. `sra` column containing the SRR accession id has to be added. In case the SRR file contains paired reads, both resulting fastq files have to be specified with the same SRR id and 'R1' and 'R2' in `read` column. SRA archive files are usually not split across lanes, so it is sufficient to specify '1' in the `lane` column. Column containing fastq files names will have to match the names of the files resulting from `fasterq-dump` tool. I haven't found a way how to specify the read suffixes, so the resulting
file name has to have a format `{srr_id}_{read}.fastq.gz` for paired reads and `{srr_id}.fastq.gz` for unpaired reads. Fastq files are automatically gzipped, using either `gzip` of multithreaded `pigz` if available on the system.

### Creating coverage plots

### Trimming

### Quality control

## Issues

