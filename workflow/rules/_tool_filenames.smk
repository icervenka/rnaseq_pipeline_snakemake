# string constants of filenames created by different programs

# Aligners ---------------------------------------------------------------------
COMMON_BAM_NAME = "aligned_sorted"

## STAR
STAR_BAM_NAME = "Aligned.sortedByCoord.out"
STAR_LOGFILES = ['Log.out', 'Log.final.out', 'Log.progress.out']

## Hisat2
HISAT_LOG_FILES = ['hisat.log']

## Tophat
TOPHAT_BAM_NAME = "accepted_hits"
TOPHAT_LOG_FILES = ["align_summary.txt"]

## Kallisto
KALLISTO_QUANT_NAME = "abundance"
KALLISTO_BAM_NAME = "pseudoalignments"
KALLISTO_LOGFILES = ["run_info.json"]

## Salmon
SALMON_QUANT_NAME = "quant"
SALMON_LOG_FILES = ['cmd_info.json', 'lib_format_counts.json']


## Samtools
SAMTOOLS_LOG_FILES = ['samtools.log']

# Read counting tools ----------------------------------------------------------
COMMON_COUNT_NAME = "count_matrix.txt"
COMMON_TRANSCRIPT_COUNT_NAME = "transcript_count_matrix.txt"

## Featurecounts
FEATURECOUNTS_COUNT_NAME = "gene_counts.txt"
FEATURECOUNTS_SUMMARY_NAME = "gene_counts.txt.summary"
FEATURECOUNTS_LOG_FILES = ["featurecounts.log"]

## HTSeq count
HTSEQ_COUNT_NAME = "gene_counts.txt"
HTSEQ_LOG_FILES = ["htseq.log"]

## Stringtie
STRINGTIE_GTF_FILE = "transcripts.gtf"
STRINGTIE_COUNT_NAME = "gene_counts.txt"
STRINGTIE_MERGED_FILE = "merged.gtf"
STRINGTIE_TPM_FILE = "stringtie_tpm.txt"
STRINGTIE_FPKM_FILE = "stringtie_fpkm.txt"
STRINGTIE_COMBINED_FILE = "stringtie_samples_combined.txt"
BALLGOWN_INPUT_FILES = ["e_data.ctab", "i_data.ctab", "t_data.ctab", "e2t.ctab", "i2t.ctab"]

## Cufflinks
CUFFLINKS_GTF_FILE = "transcripts.gtf"
CUFFLINKS_MERGED_FILE = "merged.gtf"

# Diffexp tools ----------------------------------------------------------------



# Result archive ---------------------------------------------------------------
RESULT_ARCHIVE_DIRS = [COVERAGE_OUTDIR, COUNT_OUTDIR, DIFFEXP_OUTDIR, LOG_DIR,
GSEA_INPUT_OUTDIR]