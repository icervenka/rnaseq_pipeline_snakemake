def get_align_output_files(wildcards):
    return expand(rules.align_out.output, sample=Samples)


def get_align_log_files(wildcards):
    return expand(rules.move_align_log.output, sample=Samples)

# TODO Coverage-search algorithm is turned on, making this step very slow
## Please try running TopHat again with the option (--no-coverage-search) if this step takes too much time or memory.
# TODO create transcriptome index for tophat only once and reuse it
# TODO a lot of log files are created in log directory, move them all

rule align:
    input:
        sample=get_fq,
        gtf=config["gtf"]
    output:
        bam=opj(ALIGN_OUTDIR, "{sample}", TOPHAT_BAM_FILE),
    params:
        metadata=Metadata,
        fastq_dir=FASTQ_CURRENT_DIR,
        align_outdir=opj(ALIGN_OUTDIR, "{sample}"),
        index=config["index"],
        stranded=lambda wildcards: stranded_param(wildcards, "tophat"),
        extra=has_extra_config(config["align"]["extra"], config_extra["align"])
    log:
        expand(opj(ALIGN_OUTDIR, "{{sample}}", "{log}"), log=TOPHAT_LOG_FILES),
    threads: 
        config["threads"]
    conda:
        CONDA_ALIGN_OTHER_ENV
    script:
        opj(CD2UP, WRAPPER_DIR, "tophat_wrapper.py")


rule align_out:
    input:
        rules.align.output.bam,
    output:
        opj(ALIGN_OUTDIR, "{sample}", COMMON_BAM_FILE)
    shell:
        "mv {input} {output}"


rule move_align_log:
    input:
        rules.align.log,
    output:
        expand(opj(ALIGN_LOG_OUTDIR, "{{sample}}", "{log}"), log=TOPHAT_LOG_FILES),
    shell:
        "mv {input} {output}"


include: "bam_index.smk"
