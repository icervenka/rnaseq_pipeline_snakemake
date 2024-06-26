import re


def get_multiqc_output_files(wildcards):
    return rules.multiqc.output.html


report_basename = re.sub("\s+", "", config["experiment_name"])


rule multiqc:
    input:
        get_trim_log_files,
        get_fastqc_output_files,
        get_align_log_files,
        get_count_log_files
    output:
        html=opj(MULTIQC_LOG_OUTDIR, report_basename + ".html"),
    params:
        name=report_basename,
        outdir=MULTIQC_LOG_OUTDIR
    conda:
        CONDA_SHARED_ENV
    script:
        opj(CD2UP, WRAPPER_DIR, "multiqc_wrapper.py")

        
