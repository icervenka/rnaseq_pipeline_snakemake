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
