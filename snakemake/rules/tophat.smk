# TODO replace with conda environment
rule align:
    input:
        sample = get_fq,
        index = config["index"],
        gtf = config["gtf"]
    output:
        bam = ALIGN_OUTDIR + "{sample}/" + TOPHAT_BAM_NAME + ".bam",
        log = ALIGN_OUTDIR + "{sample}/align_summary.txt"
    params:
        extra = config['align']['tophat_extra']
    threads:
        config["threads"]
    run:
        shell(
            "export PS1=; "
            "source /usr/local/bin/miniconda3/etc/profile.d/conda.sh; "
            "conda activate tophat; "
            "tophat2 "
            "{params.extra} "
            "-p {threads} "
            "-G {input.gtf} "
            "-o {ALIGN_OUTDIR} "
            "{input.index} "
            "{input.sample} "
        )

rule rename_bam:
    input:
        rules.align.output.bam
    output:
        ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".bam"
    shell:
        "mv {input} {output}"

rule move_log:
    input:
        rules.align.output.log
    output:
        ALIGN_LOG_OUTDIR + "{sample}/align_summary.txt"
    shell:
        "mv {input} {output}"

rule all_align:
    input:
        expand(ALIGN_OUTDIR + "{sample}/" +
               COMMON_BAM_NAME + ".bam", sample=Samples),
        expand(ALIGN_LOG_OUTDIR + "{sample}/align_summary.txt", sample=Samples)
    output:
        touch(LOG_DIR + "align.completed")
