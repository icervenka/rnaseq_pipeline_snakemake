def get_align_output_files(wildcards):
    return expand(ALIGN_OUTDIR+"{sample}/"+COMMON_BAM_NAME+".bam", sample=Samples)

def get_align_log_files(wildcards):
    return expand(ALIGN_LOG_OUTDIR+"{sample}/{log}", sample=Samples, log=STAR_LOGFILES)

rule align:
    input:
        sample=get_fq,
        index=config["index"]
    output:
        bam=temp(ALIGN_OUTDIR + "{sample}/" + STAR_BAM_NAME + ".bam"),
        log=expand(ALIGN_OUTDIR + "{{sample}}/{log}", log=STAR_LOGFILES),
    params:
        extra=config_extra['align'][config['align']['extra']],
        metadata=Metadata,
        fastq_dir=FASTQ_INPUT_DIR
    threads:
        config["threads"]
    conda:
        CONDA_ALIGN_GENERAL_ENV
    script:
        "../scripts/star_wrapper.py"

rule rename_bam:
    input:
        rules.align.output.bam
    output:
        ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".bam"
    shell:
        "mv {input} {output}"

rule move_align_log:
    input:
        rules.align.output.log
    output:
        expand(ALIGN_LOG_OUTDIR + "{{sample}}/{log}", log=STAR_LOGFILES)
    params:
        outdir=ALIGN_LOG_OUTDIR + "{sample}"
    shell:
        "mv {input} {params.outdir}"

include: "bam_index.smk"
