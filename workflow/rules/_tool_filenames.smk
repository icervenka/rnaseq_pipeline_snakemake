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