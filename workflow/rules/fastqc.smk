def get_fastqc_output_files(wildcards):
    return(
        expand(rules.fastqc.output.html, sample=Metadata.fq) + 
        expand(rules.fastqc.output.zip, sample=Metadata.fq)
    ) + get_fastqc_processed_output_files(wildcards)


rule fastqc:
    input:
        opj(FASTQ_DIR, "{sample}")
    output:
        html=opj(FASTQC_LOG_OUTDIR, "{sample}.html"),
        # the suffix _fastqc.zip is necessary for multiqc to find the file.
        # If not using multiqc, you are free to choose an arbitrary filename
        zip=opj(FASTQC_LOG_OUTDIR, "{sample}_fastqc.zip")
    params:
        ""
    threads:
        1
    log:
        opj(FASTQC_LOG_OUTDIR, "{sample}.log")
    conda:
        CONDA_SHARED_ENV
    wrapper:
        "v1.21.0/bio/fastqc"

if FASTQ_DIR != FASTQ_CURRENT_DIR:

    def get_fastqc_processed_output_files(wildcards):
        return(
            expand(rules.fastqc_processed.output.html, sample=Metadata.fq) + 
            expand(rules.fastqc_processed.output.zip, sample=Metadata.fq)
        )

    rule fastqc_processed:
        input:
            opj(FASTQ_CURRENT_DIR, "{sample}")
        output:
            html=opj(FASTQ_CURRENT_DIR, "{sample}.html"),
            # the suffix _fastqc.zip is necessary for multiqc to find the file.
            # If not using multiqc, you are free to choose an arbitrary filename
            zip=opj(FASTQC_LOG_OUTDIR, "processed", "{sample}_fastqc.zip")
        params:
            ""
        threads:
            1
        log:
            opj(FASTQC_LOG_OUTDIR, "processed", "{sample}.log")
        conda:
            CONDA_SHARED_ENV
        wrapper:
            "v1.21.0/bio/fastqc"

else: 
    def get_fastqc_processed_output_files(wildcards):
        return []
