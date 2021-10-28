# folder organization of the pipelline input and output files
RULES_DIR = "snakemake/rules/"
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

# string constants of filenames created by different programs
STAR_BAM_NAME = "Aligned.sortedByCoord.out"
STAR_LOGFILES = ['Log.out', 'Log.final.out', 'Log.progress.out']
TOPHAT_BAM_NAME = "accepted_hits"
KALLISTO_BAM_NAME = "abundance"
COMMON_BAM_NAME = "aligned_sorted"
BALLGOWN_INPUT_FILES = ["e_data", "i_data", "t_data", "e2t", "i2t"]
