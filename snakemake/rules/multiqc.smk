import re

rule multiqc:
    input:
        expand(LOG_DIR + "qc/{file}.html", file=Metadata.fq),
        expand(LOG_DIR + "qc/{file}_fastqc.zip", file=Metadata.fq),
        get_align_log_files,
        get_count_log_files,
    output:
        html=LOG_DIR + "qc/" + re.sub("\s+", "", config["experiment_name"]) + ".html",
    params:
        name=re.sub("\s+", "", config["experiment_name"]),
        outdir=LOG_DIR + "qc/"
    shell:
        "multiqc -f -o {params.outdir} -n {params.name} {params.outdir}"
