# folder organization of the pipelline input and output files
RULES_DIR = "workflow/rules/"
ENV_DIR = "envs/"
SRA_DIR = "sra/"
FASTQ_DIR = "fastq/"
TRIMMED_DIR = "trimmed/"
LOG_DIR = "logs/"
ALIGN_OUTDIR = "align/"
ALIGN_LOG_OUTDIR= LOG_DIR + ALIGN_OUTDIR
COVERAGE_OUTDIR = "coverage/"
COUNT_OUTDIR = "counts/"
MERGE_OUTDIR = "merged/"
COUNT_LOG_OUTDIR = LOG_DIR + COUNT_OUTDIR
DIFFEXP_OUTDIR = "diffexp/"
DIFFEXP_LOG_OUTDIR = LOG_DIR + DIFFEXP_OUTDIR
GSEA_INPUT_OUTDIR = "gsea_input/"
ARCHIVE_DIR = "archive/"

# string constants of filenames created by different programs
STAR_BAM_NAME = "Aligned.sortedByCoord.out"
STAR_LOGFILES = ['Log.out', 'Log.final.out', 'Log.progress.out']
TOPHAT_BAM_NAME = "accepted_hits"
KALLISTO_QUANT_NAME = "abundance"
KALLISTO_BAM_NAME = "pseudocounts"
SALMON_QUANT_NAME = "quant"
SALMON_LOG_FILES = ['cmd_info.json', 'lib_format_counts.json']
COMMON_BAM_NAME = "aligned_sorted"
BALLGOWN_INPUT_FILES = ["e_data", "i_data", "t_data", "e2t", "i2t"]
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
CONDA_DIFFEXP_GENERAL_ENV = ENV_DIR + 'diffexp_general.yaml'