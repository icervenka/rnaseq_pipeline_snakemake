def get_align_output_files(wildcards):
    return (
        expand(rules.align_out.output, sample=Samples) +
        expand(rules.align.output.splicesite, sample=Samples)
    )


def get_align_log_files(wildcards):
    return (
        expand(rules.align.log.log, sample=Samples) + 
        expand(rules.align.log.met, sample=Samples)
    )


# TODO rule for generating splice sites


rule align:
    input:
        sample=get_fq
    output:
        sam=temp(opj(ALIGN_OUTDIR, "{sample}", COMMON_BAM_NAME + ".sam")),
        splicesite=opj(ALIGN_OUTDIR, "{sample}", "splice_sites.txt"),  
    params:
        metadata=Metadata,
        fastq_dir=FASTQ_INPUT_DIR,
        index=config["index"],
        extra=has_extra_config(config["align"]["extra"], config_extra["align"]),
    log:
        log=expand(opj(ALIGN_LOG_OUTDIR, "{{sample}}", "{log}"), log=HISAT_LOG_FILES),
        met=expand(opj(ALIGN_LOG_OUTDIR, "{{sample}}", "{log}"), log=["hisat_metrics.log"]),
    threads: 
        config["threads"]
    conda:
        CONDA_ALIGN_GENERAL_ENV
    script:
        "../scripts/hisat_wrapper.py"


rule align_out:
    input:
        rules.align.output.sam,
    output:
        opj(ALIGN_OUTDIR, "{sample}", COMMON_BAM_NAME + ".bam"),
    params:
        compression=9,
    threads: config["threads"]
    conda:
        CONDA_ALIGN_GENERAL_ENV
    shell:
        """
        samtools view  \
        -uSh \
        {input} | \
        samtools sort \
        -l {params.compression}  \
        -@ {threads} \
        --output-fmt BAM \
        -o {output}
        """


include: "bam_index.smk"
