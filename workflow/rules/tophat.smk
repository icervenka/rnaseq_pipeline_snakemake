def get_align_output_files(wildcards):
    return expand(ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".bam", sample=Samples)


def get_align_log_files(wildcards):
    return expand(ALIGN_LOG_OUTDIR + "{sample}/align_summary.txt", sample=Samples)

# TODO add library type
# TODO function to create align_outdir in params
# TODO Coverage-search algorithm is turned on, making this step very slow
## Please try running TopHat again with the option (--no-coverage-search) if this step takes too much time or memory.
# TODO create transcriptome index for tophat only once and reuse it
# TODO a lot of log files are created in log directory, move them all

rule align:
    input:
        sample=get_fq,
        gtf=config["gtf"]
    output:
        bam=ALIGN_OUTDIR + "{sample}/" + TOPHAT_BAM_NAME + ".bam",
        log=ALIGN_OUTDIR + "{sample}/align_summary.txt",
    params:
        metadata=Metadata,
        fastq_dir=FASTQ_INPUT_DIR,
        align_outdir=ALIGN_OUTDIR + "{sample}/",
        index=config["index"],
        extra=has_extra_config(config["align"]["extra"], config_extra["align"])
    threads: config["threads"]
    conda:
        CONDA_ALIGN_OTHER_ENV
    script:
        "../scripts/tophat_wrapper.py"


rule align_out:
    input:
        rules.align.output.bam,
    output:
        ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".bam",
    shell:
        "mv {input} {output}"


rule move_log:
    input:
        rules.align.output.log,
    output:
        ALIGN_LOG_OUTDIR + "{sample}/align_summary.txt",
    shell:
        "mv {input} {output}"


include: "bam_index.smk"
