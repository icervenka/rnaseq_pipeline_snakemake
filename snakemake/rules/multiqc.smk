import re

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
