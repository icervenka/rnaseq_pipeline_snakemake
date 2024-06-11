# folder organization of the pipelline input and output files

_ds_new = {
    "INPUT_DIR": "input",
    "OUTPUT_DIR": "output",
    "WORKFLOW_DIR": "workflow"
}

_ds = {
    "CD1UP": "../",
    "CD2UP": "../../",
    "CD3UP": "../../../",
    "ENV_DIR": "workflow/envs/",
    "SCRIPT_DIR": "workflow/scripts/",
    "SRA_DIR": "sra/",
    "FASTQ_DIR": "fastq/",
    "FASTQ_CURRENT_DIR": "fastq/",
    "FASTQ_PREPROCESSED_DIR": "fastq_preprocessed/",
    "FASTQ_TRIMMED_DIR": "fastq_trimmed/",
    "ALIGN_OUTDIR": "align/",
    "COVERAGE_OUTDIR": "coverage/",
    "COUNT_OUTDIR": "counts/",
    "DIFFEXP_OUTDIR": "diffexp/",
    "GSEA_INPUT_OUTDIR": "misc/gsea_input/",
    "QC_DIR": "qc",
    "LOG_DIR": "logs/",
    "ARCHIVE_DIR": "archive/"
}

# folder organization of the pipelline log_files
_ds.update({
    "SUBSAMPLE_LOG_OUTDIR": opj(_ds["LOG_DIR"], "subsample/"),
    "TRIM_LOG_OUTDIR": opj(_ds["LOG_DIR"], "trim/"),
    "ALIGN_LOG_OUTDIR": opj(_ds["LOG_DIR"], "align/"),
    "COUNT_LOG_OUTDIR": opj(_ds["LOG_DIR"], "counts/"),
    "DIFFEXP_LOG_OUTDIR": opj(_ds["LOG_DIR"], "diffexp/"),
    "FASTQC_LOG_OUTDIR": opj(_ds["LOG_DIR"], "fastqc/"),
    "MULTIQC_LOG_OUTDIR":opj( _ds["LOG_DIR"], "multiqc/")
})

# folder organization of other folders that depend on the default organization
_ds.update({
    "DIFFEXP_ANALYSIS_OUTDIR": opj(_ds["DIFFEXP_OUTDIR"], config["diffexp"]["outdir"]),
    "CUFFCOMPARE_OUTDIR": opj(_ds["COUNT_OUTDIR"], "cuffcompare"),
    "CUFFNORM_OUTDIR": opj(_ds["COUNT_OUTDIR"], "cuffnorm")

})

# folder organization of the diffexp files
_ds.update({
    "DEGFILES_OUTDIR": opj(_ds["DIFFEXP_ANALYSIS_OUTDIR"], "degfiles/"),
    "REPORTS_OUTDIR": opj(_ds["DIFFEXP_ANALYSIS_OUTDIR"], "reports/"),
    "RDS_OUTDIR": opj(_ds["DIFFEXP_ANALYSIS_OUTDIR"], "rds/"),
    "ANALYSIS_PARAM_OUTDIR": opj(_ds["DIFFEXP_ANALYSIS_OUTDIR"], "analysis_params/")
})

# conda environments
# conda environments are invoked from rules dir, and directory structure is organized
# in relation to base Snakemake dir, going up 2 dirs has to be prepended to the path
_envs = {
    "CONDA_SHARED_ENV":  opj(_ds["CD2UP"], _ds["ENV_DIR"], "shared.yaml"),
    "CONDA_SRA_TOOLS_ENV":  opj(_ds["CD2UP"], _ds["ENV_DIR"], "sra_tools.yaml"),
    "CONDA_PREPROCESS_ENV":  opj(_ds["CD2UP"], _ds["ENV_DIR"], "preprocess.yaml"),
    "CONDA_CUTADAPT_ENV":  opj(_ds["CD2UP"], _ds["ENV_DIR"], "cutadapt.yaml"),
    "CONDA_QC_ENV":  opj(_ds["CD2UP"], _ds["ENV_DIR"], "qc.yaml"),
    "CONDA_ALIGN_GENERAL_ENV":  opj(_ds["CD2UP"], _ds["ENV_DIR"], "align_general.yaml"),
    "CONDA_ALIGN_OTHER_ENV": opj(_ds["CD2UP"], _ds["ENV_DIR"], "align_other.yaml"),
    "CONDA_COUNT_GENERAL_ENV": opj(_ds["CD2UP"], _ds["ENV_DIR"], "count_general.yaml"),
    "CONDA_COUNT_CUFFLINKS_ENV":  opj(_ds["CD2UP"], _ds["ENV_DIR"], "count_cufflinks.yaml"),
    "CONDA_R_GENERAL_ENV":  opj(_ds["CD2UP"], _ds["ENV_DIR"], "r_general.yaml"),
    "CONDA_DIFFEXP_GENERAL_ENV":  opj(_ds["CD2UP"], _ds["ENV_DIR"], "diffexp_general.yaml")
}

# TODO implement
# # Result archive ---------------------------------------------------------------
# RESULT_ARCHIVE_DIRS = [COVERAGE_OUTDIR, COUNT_OUTDIR, DIFFEXP_OUTDIR, LOG_DIR,
# GSEA_INPUT_OUTDIR]

# update into global environment
globals().update(_ds)
globals().update(_envs)