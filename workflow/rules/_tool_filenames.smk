# string constants of filenames created by different programs
# NOTES
# - files are included with extension
# - format is roughly <TOOL>_<TYPE>_<SUBTYPES>_FILES(S)
# - prefix COMMON is used to unify filenames between tools that are then 
#   consumed by other steps
# - ending with _FILE denotes a string, ending with _FILES denotes array of 
#   strings
# - full cuffcompare files are specified in the command with an outprefix

_tf = {
    
    # Aligners ---------------------------------------------------------------------
    "COMMON_BAM_FILE": "aligned_sorted.bam",

    ## STAR
    "STAR_BAM_FILE": "Aligned.sortedByCoord.out.bam",
    "STAR_LOGFILES": ['Log.out', 'Log.final.out', 'Log.progress.out'],

    ## Hisat2
    "HISAT_BAM_FILE": "aligned_sorted.bam",
    "HISAT_SAM_FILE": "aligned_sorted.sam",
    "HISAT_LOG_FILES": ['hisat.log'],

    ## Tophat
    "TOPHAT_BAM_FILE": "accepted_hits.bam",
    "TOPHAT_LOG_FILES": ["align_summary.txt"],

    ## Kallisto
    "KALLISTO_QUANT_H5_FILE": "abundance.h5",
    "KALLISTO_QUANT_TSV_FILE": "abundance.tsv",
    "KALLISTO_BAM_FILE": "pseudoalignments.bam",
    "KALLISTO_LOGFILES": ["run_info.json"],

    ## Salmon
    "SALMON_QUANT_FILE": "quant.sf",
    "SALMON_CMDLOG_FILES": ['salmon.log'],
    "SALMON_LOG_FILES": ['cmd_info.json', 'lib_format_counts.json'],


    ## Samtools
    "SAMTOOLS_LOG_FILES": ['samtools.log'],

    # Read counting tools ----------------------------------------------------------
    "COMMON_COUNT_FILE": "count_matrix.txt",
    "COMMON_TRANSCRIPT_COUNT_FILE": "transcript_count_matrix.txt",

    ## Featurecounts
    "FEATURECOUNTS_COUNT_FILE": "gene_counts.txt",
    "FEATURECOUNTS_SUMMARY_FILE": "gene_counts.txt.summary",
    "FEATURECOUNTS_LOG_FILES": ["featurecounts.log"],

    ## HTSeq count
    "HTSEQ_COUNT_FILE": "gene_counts.txt",
    "HTSEQ_LOG_FILES": ["htseq.log"],

    ## Stringtie
    "STRINGTIE_GTF_FILE": "transcripts.gtf",
    "STRINGTIE_COUNT_FILE": "gene_counts.txt",
    "STRINGTIE_MERGED_FILE": "merged.gtf",
    "STRINGTIE_TPM_FILE": "stringtie_tpm.txt",
    "STRINGTIE_FPKM_FILE": "stringtie_fpkm.txt",
    "STRINGTIE_COMBINED_FILE": "stringtie_samples_combined.txt",
    "BALLGOWN_INPUT_FILES": ["e_data.ctab", "i_data.ctab", "t_data.ctab", "e2t.ctab", "i2t.ctab"],
    "STRINGTIE_LOG_FILES": ["stringtie.log"],
    "STRINGTIE_MERGE_LOG_FILES": ["stringtie_merge.log"],
    "STRINGTIE_COUNT_LOG_FILES": ["stringtie_quant.log"],

    ## Cufflinks
    "CUFFLINKS_GTF_FILE": "transcripts.gtf",
    "CUFFLINKS_MERGED_FILE": "merged.gtf",
    "CUFFCOMPARE_FILES": [
        ".loci", ".stats", ".tracking"
    ],
    "CUFFCOMPARE_GTF_FILE": ".combined.gtf",
    "CUFFNORM_COUNT_TABLE_FILE": "genes.count_table",
    "CUFFNORM_ATTR_TABLE_FILE": "genes.attr_table",
    "CUFFNORM_SAMPLE_TABLE_FILE": "samples.table",
    "CUFFNORM_COUNT_FILES":  [
        "cds.attr_table", "cds.count_table", "cds.fpkm_table", "genes.attr_table",  
        "genes.count_table", "genes.fpkm_table", "isoforms.attr_table", 
        "isoforms.count_table", "isoforms.fpkm_table", "samples.table", 
        "tss_groups.attr_table", "tss_groups.count_table", "tss_groups.fpkm_table"
    ],
    "CUFFLINKS_LOG_FILES": ["cufflinks.log"],
    "CUFFMERGE_RUN_LOG_FILES": ["cuffmerge_run.log"],
    "CUFFMERGE_LOG_FILES": ["cuffmerge.log"],
    "CUFFQUANT_LOG_FILES": ["cuffquant.log"],
    "CUFFNORM_LOG_FILES": ["cuffnorm.log"],

    # Diffexp tools ------------------------------------------------------------
    "CUFFDIFF_GENE_DEG_FILE": "gene_exp.diff",
    "CUFFDIFF_DIFFEXP_FILES": [
        "bias_params.info", "cds.count_tracking", "cds.diff", "cds.fpkm_tracking", 
        "cds.read_group_tracking", "cds_exp.diff", "gene_exp.diff", 
        "genes.count_tracking", "genes.fpkm_tracking", "genes.read_group_tracking", 
        "isoform_exp.diff", "isoforms.count_tracking", "isoforms.fpkm_tracking", 
        "isoforms.read_group_tracking", "promoters.diff", "read_groups.info", 
        "run.info", "splicing.diff", "tss_group_exp.diff", "tss_groups.count_tracking", 
        "tss_groups.fpkm_tracking", "tss_groups.read_group_tracking", "var_model.info"
    ],
    "CUFFDIFF_LOG_FILES": ["cuffdiff.log"]

}

globals().update(_tf)