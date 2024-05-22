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
        sam=temp(ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".sam"),
        splicesite=ALIGN_OUTDIR + "{sample}/" + "splice_sites.txt",
        
    params:
        metadata=Metadata,
        index=config["index"],
        fastq_dir=FASTQ_INPUT_DIR,
        extra=has_extra_config(config["align"]["extra"], config_extra["align"]),
    log:
        log=ALIGN_LOG_OUTDIR + "{sample}/hisat.log",
        met=ALIGN_LOG_OUTDIR + "{sample}/hisat_metrics.txt"
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
