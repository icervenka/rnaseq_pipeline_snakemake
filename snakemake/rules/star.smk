rule align:
    input:
        sample=get_fq,
        index=config["index"]
    output:
        bam=temp(ALIGN_OUTDIR + "{sample}/" + STAR_BAM_NAME + ".bam"),
        log=expand(ALIGN_OUTDIR + "{{sample}}/{log}", log=STAR_LOGFILES),
    params:
        extra=extra_config['align']['star_extra'],
        metadata=Metadata,
        fastq_dir = FASTQ_DIR
    threads:
        config["threads"]
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

rule all_align:
    input:
        expand(ALIGN_OUTDIR+"{sample}/"+COMMON_BAM_NAME+".bam", sample=Samples),
        expand(ALIGN_LOG_OUTDIR+"{sample}/{log}", sample=Samples, log=STAR_LOGFILES)
    output:
        touch(LOG_DIR + "align.completed")
