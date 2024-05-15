import re

def get_multiqc_output_files(wildcards):
    return(LOG_DIR +
           "multiqc/" +
           re.sub("\s+", "", config["experiment_name"]) +
           ".html")

rule multiqc:
    input:
        get_trim_log_files,
        get_fastqc_output_files,
        get_align_log_files,
        get_count_log_files
    output:
        html = LOG_DIR + "multiqc/" + \
            re.sub("\s+", "", config["experiment_name"]) + ".html",
    params:
        name = re.sub("\s+", "", config["experiment_name"]),
        outdir = LOG_DIR + "multiqc/"
    # TODO conda directive doesn't work with run directive
    conda:
        CONDA_SHARED_ENV
    script:
        "../scripts/multiqc_wrapper.py"

        
