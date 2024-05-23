# TODO think about refactoring with os.path.join into dictionaries or maybe yaml files
# TODO separate into folder structure and other global variables

# folder organization of the pipelline input and output files
RULES_DIR = "workflow/rules/"
ENV_DIR = "envs/"
SRA_DIR = "sra/"
FASTQ_DIR = "fastq/"
FASTQ_ALIGN_DIR = "fastq/"
FASTQ_PREPROCESSED_DIR = "fastq_preprocessed/"
FASTQ_TRIMMED_DIR = "fastq_trimmed/"
LOG_DIR = "logs/"
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

# string constants of filenames created by different programs
STAR_BAM_NAME = "Aligned.sortedByCoord.out"
STAR_LOGFILES = ['Log.out', 'Log.final.out', 'Log.progress.out']
TOPHAT_BAM_NAME = "accepted_hits"
TOPHAT_LOG_FILES = ["align_summary.txt"]
KALLISTO_QUANT_NAME = "abundance"
KALLISTO_BAM_NAME = "pseudoalignments"
KALLISTO_LOGFILES = ["run_info.json"]
SALMON_QUANT_NAME = "quant"
SALMON_LOG_FILES = ['cmd_info.json', 'lib_format_counts.json']
HISAT_LOG_FILES = ['hisat.log']

COMMON_BAM_NAME = "aligned_sorted"
BALLGOWN_INPUT_FILES = ["e_data", "i_data", "t_data", "e2t", "i2t"]

COMMON_COUNT_NAME = "counts.tsv"
FEATURECOUNTS_COUNT_NAME = "gene_counts.txt"
FEATURECOUNTS_SUMMARY_NAME = "gene_counts.txt.summary"
FEATURECOUNTS_LOG_FILES = ["featurecounts.log"]
HTSEQ_COUNT_NAME = "gene_counts.txt"
HTSEQ_LOG_FILES = ["htseq.log"]

RESULT_ARCHIVE_DIRS = [COVERAGE_OUTDIR, COUNT_OUTDIR, DIFFEXP_OUTDIR, LOG_DIR,
GSEA_INPUT_OUTDIR]

# Conda environments
CONDA_SHARED_ENV = ENV_DIR + 'shared.yaml'
CONDA_SRA_TOOLS_ENV = ENV_DIR + 'sra_tools.yaml'
CONDA_PREPROCESS_ENV = ENV_DIR + 'preprocess.yaml'
CONDA_CUTADAPT_ENV = ENV_DIR + 'cutadapt.yaml'
CONDA_QC_ENV = ENV_DIR + 'qc.yaml'
CONDA_ALIGN_GENERAL_ENV = ENV_DIR + 'align_general.yaml'
CONDA_ALIGN_OTHER_ENV = ENV_DIR + 'align_other.yaml'
CONDA_COUNT_GENERAL_ENV = ENV_DIR + 'count_general.yaml'
CONDA_COUNT_CUFFLINKS_ENV = ENV_DIR + 'count_cufflinks.yaml'
CONDA_R_GENERAL_ENV = ENV_DIR + 'r_general.yaml'
CONDA_DIFFEXP_GENERAL_ENV = ENV_DIR + 'diffexp_general.yaml'