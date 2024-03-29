def get_fastqc_output_files(wildcards):
    return(
        expand(LOG_DIR + "fastqc/{sample}.html", sample=Metadata.fq) + 
        expand(LOG_DIR + "fastqc/{sample}_fastqc.zip", sample=Metadata.fq)
    )

rule fastqc:
    input:
        FASTQ_DIR + "{sample}"
    output:
        html=LOG_DIR + "fastqc/{sample}.html",
        # the suffix _fastqc.zip is necessary for multiqc to find the file.
        # If not using multiqc, you are free to choose an arbitrary filename
        zip=LOG_DIR + "fastqc/{sample}_fastqc.zip"
    params:
        ""
    threads:
        1
    log:
        LOG_DIR + "fastqc/{sample}.log"
    wrapper:
        "v1.21.0/bio/fastqc"
