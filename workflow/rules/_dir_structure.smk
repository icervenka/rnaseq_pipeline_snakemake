# TODO think about refactoring with os.path.join into dictionaries or maybe yaml files
# TODO separate into folder structure and other global variables

# folder organization of the pipelline input and output files
RULES_DIR = "workflow/rules/"
ENV_DIR = "envs/"
SRA_DIR = "sra/"
FASTQ_DIR = "fastq/"
FASTQ_ALIGN_DIR = "fastq/"
FASTQ_CURRENT_DIR = "fastq/"
FASTQ_PREPROCESSED_DIR = "fastq_preprocessed/"
FASTQ_TRIMMED_DIR = "fastq_trimmed/"
LOG_DIR = "logs/"
SUBSAMPLE_LOG_OUTDIR = LOG_DIR + "subsample/"
TRIM_LOG_OUTDIR = LOG_DIR + "trim/"
ALIGN_OUTDIR = "align/"
ALIGN_LOG_OUTDIR= LOG_DIR + ALIGN_OUTDIR
COVERAGE_OUTDIR = "coverage/"
COUNT_OUTDIR = "counts/"
MERGE_OUTDIR = "merged/"
COUNT_LOG_OUTDIR = LOG_DIR + COUNT_OUTDIR
DIFFEXP_OUTDIR = "diffexp/"
DIFFEXP_LOG_OUTDIR = LOG_DIR + DIFFEXP_OUTDIR
GSEA_INPUT_OUTDIR = "misc/gsea_input/"
ARCHIVE_DIR = "archive/"