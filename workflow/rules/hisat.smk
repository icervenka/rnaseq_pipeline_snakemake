def get_align_output_files(wildcards):
    return expand(
        ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".bam", sample=Samples
    ) + expand(ALIGN_OUTDIR + "{sample}/" + "splice_sites.txt", sample=Samples)


def get_align_log_files(wildcards):
    return expand(ALIGN_LOG_OUTDIR + "{sample}/hisat.log", sample=Samples) + expand(
        ALIGN_LOG_OUTDIR + "{sample}/hisat_metrics.txt", sample=Samples
    )


# TODO rule for generating splice sites


rule align:
    input:
        sample=get_fq,
    output:
        sam=temp(ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".sam"),
        splicesite=ALIGN_OUTDIR + "{sample}/" + "splice_sites.txt",
        log=ALIGN_LOG_OUTDIR + "{sample}/hisat.log",
        met=ALIGN_LOG_OUTDIR + "{sample}/hisat_metrics.txt",
    params:
        metadata=Metadata,
        index=config["index"],
        fastq_dir=FASTQ_INPUT_DIR,
        extra=has_extra_config(config["align"]["extra"], config_extra["align"]),
    threads: config["threads"]
    conda:
        CONDA_ALIGN_GENERAL_ENV
    script:
        "../scripts/hisat_wrapper.py"


rule rename_bam:
    input:
        rules.align.output.sam,
    output:
        ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".bam",
    params:
        compression=9,
        flag="0x2",
    threads: config["threads"]
    conda:
        CONDA_ALIGN_GENERAL_ENV
    shell:
        """
        samtools view  -uSh -f {params.flag}  {input} | \
        samtools sort -l {params.compression}  -@ {threads} \
        --output-fmt BAM -o {output}
        """


include: "bam_index.smk"
