def get_align_output_files(wildcards):
    return expand(rules.align_out.output, sample=Samples)

def get_align_log_files(wildcards):
    return expand(rules.move_align_log.output, sample=Samples)

rule align:
    input:
        sample=get_fq
    output:
        bam=opj(ALIGN_OUTDIR, "{sample}", STAR_BAM_FILE)
    params:
        metadata=Metadata,
        fastq_dir=FASTQ_CURRENT_DIR,
        index=config["index"],
        extra=has_extra_config(config["align"]["extra"], config_extra["align"])
    log:
        expand(opj(ALIGN_OUTDIR, "{{sample}}", "{log}"), log=STAR_LOGFILES)
    threads:
        config["threads"]
    conda:
        CONDA_ALIGN_GENERAL_ENV
    script:
        opj(CD2UP, WRAPPER_DIR, "star_wrapper.py")


rule align_out:
    input:
        rules.align.output.bam
    output:
        opj(ALIGN_OUTDIR, "{sample}", COMMON_BAM_FILE)
    shell:
        "mv {input} {output}"

rule move_align_log:
    input:
        rules.align.log
    output:
        expand(opj(ALIGN_LOG_OUTDIR, "{{sample}}", "{log}"), log=STAR_LOGFILES)
    params:
        outdir=opj(ALIGN_LOG_OUTDIR, "{sample}")
    shell:
        "mv {input} {params.outdir}"

include: "bam_index.smk"
