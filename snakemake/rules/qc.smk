import re

rule fastqc:
    input:
        lambda wildcards: Metadata.fq[Metadata.fq == wildcards.fq]
    output:
        html=LOG_DIR + "qc/{fq}.html",
        # the suffix _fastqc.zip is necessary for multiqc to find the file.
        # If not using multiqc, you are free to choose an arbitrary filename
        zip=LOG_DIR + "qc/{fq}_fastqc.zip"
    params:
        ""
    threads:
        config['threads']
    log:
        LOG_DIR + "qc/{fq}.log"
    wrapper:
        "0.40.0/bio/fastqc"

rule multiqc:
    input:
        expand(rules.fastqc.output.html, fq=Fqs),
        expand(rules.fastqc.output.zip, fq=Fqs),
        expand(rules.move_align_log.output[0], sample=Samples),
        expand(rules.move_count_log.output[0], sample=Samples)
    output:
        html=LOG_DIR + "qc/" + re.sub("\s+", "", config["experiment_name"]) + ".html",
    params:
        name=re.sub("\s+", "", config["experiment_name"]),
        outdir=LOG_DIR + "qc/"
    shell:
        "multiqc -f -o {params.outdir} -n {params.name} {params.outdir}"

rule qc_all:
    input:
        expand(rules.fastqc.output.html, fq=Fqs),
        expand(rules.fastqc.output.zip, fq=Fqs),
        rules.multiqc.output.html
    output:
        touch(LOG_DIR + "qc.completed")
