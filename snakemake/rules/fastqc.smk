rule fastqc:
    input:
        get_fq
    output:
        html=LOG_DIR + "qc/{sample}.html",
        # the suffix _fastqc.zip is necessary for multiqc to find the file.
        # If not using multiqc, you are free to choose an arbitrary filename
        zip=LOG_DIR + "qc/{sample}_fastqc.zip"
    params:
        ""
    threads:
        config['threads']
    log:
        LOG_DIR + "qc/{sample}.log"
    wrapper:
        "0.40.0/bio/fastqc"
