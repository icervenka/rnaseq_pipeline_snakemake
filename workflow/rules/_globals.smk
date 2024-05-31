# folder organization of the pipelline input and output files
_ds = {
    "RULES_DIR": "workflow/rules/",
    "SCRIPT_DIR": "workflow/scripts/",
    "ENV_DIR": "envs/",
    "SRA_DIR": "sra/",
    "FASTQ_DIR": "fastq/",
    "FASTQ_CURRENT_DIR": "fastq/",
    "FASTQ_PREPROCESSED_DIR": "fastq_preprocessed/",
    "FASTQ_TRIMMED_DIR": "fastq_trimmed/",
    "ALIGN_OUTDIR": "align/",
    "COVERAGE_OUTDIR": "coverage/",
    "COUNT_OUTDIR": "counts/",
    "MERGE_OUTDIR": "merged/",
    "DIFFEXP_OUTDIR": "diffexp/",
    "GSEA_INPUT_OUTDIR": "misc/gsea_input/",
    "LOG_DIR": "logs/",
    "ARCHIVE_DIR": "archive/",
    "CD1UP": "../",
    "CD2UP": "../../",
    "CD3UP": "../../../"
}

# folder organization of the pipelline log_files
_ds_log = {
    "SUBSAMPLE_LOG_OUTDIR": opj(_ds["LOG_DIR"], "subsample/"),
    "TRIM_LOG_OUTDIR": opj(_ds["LOG_DIR"], "trim/"),
    "ALIGN_LOG_OUTDIR": opj(_ds["LOG_DIR"], "align/"),
    "COUNT_LOG_OUTDIR": opj(_ds["LOG_DIR"], "counts/"),
    "DIFFEXP_LOG_OUTDIR": opj(_ds["LOG_DIR"], "diffexp/"),
    "FASTQC_LOG_OUTDIR": opj(_ds["LOG_DIR"], "fastqc/"),
    "MULTIQC_LOG_OUTDIR":opj( _ds["LOG_DIR"], "multiqc/")
}

# conda environments
_envs = {
    "CONDA_SHARED_ENV":  opj(_ds["ENV_DIR"], "shared.yaml"),
    "CONDA_SRA_TOOLS_ENV":  opj(_ds["ENV_DIR"], "sra_tools.yaml"),
    "CONDA_PREPROCESS_ENV":  opj(_ds["ENV_DIR"], "preprocess.yaml"),
    "CONDA_CUTADAPT_ENV":  opj(_ds["ENV_DIR"], "cutadapt.yaml"),
    "CONDA_QC_ENV":  opj(_ds["ENV_DIR"], "qc.yaml"),
    "CONDA_ALIGN_GENERAL_ENV":  opj(_ds["ENV_DIR"], "align_general.yaml"),
    "CONDA_ALIGN_OTHER_ENV": opj(_ds["ENV_DIR"], "align_other.yaml"),
    "CONDA_COUNT_GENERAL_ENV": opj( _ds["ENV_DIR"], "count_general.yaml"),
    "CONDA_COUNT_CUFFLINKS_ENV":  opj(_ds["ENV_DIR"], "count_cufflinks.yaml"),
    "CONDA_R_GENERAL_ENV":  opj(_ds["ENV_DIR"], "r_general.yaml"),
    "CONDA_DIFFEXP_GENERAL_ENV":  opj(_ds["ENV_DIR"], "diffexp_general.yaml")
}

# update into global environment
globals().update(_ds)
globals().update(_ds_log)
globals().update(_envs)