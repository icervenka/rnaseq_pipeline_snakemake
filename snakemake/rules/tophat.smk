def get_align_output_files(wildcards):
    return expand(ALIGN_OUTDIR + "{sample}/" +
        COMMON_BAM_NAME + ".bam", sample=Samples),

def get_align_log_files(wildcards):
    return expand(ALIGN_LOG_OUTDIR + "{sample}/align_summary.txt", sample=Samples)

# TODO replace with conda environment
rule align:
    input:
        sample=get_fq,
        index=config["index"],
        gtf=config["gtf"]
    output:
        bam=ALIGN_OUTDIR + "{sample}/" + TOPHAT_BAM_NAME + ".bam",
        log=ALIGN_OUTDIR + "{sample}/align_summary.txt"
    params:
        extra=config_extra['align']['tophat_extra'],
        metadata=Metadata,
        fastq_dir=FASTQ_DIR
    threads:
        config["threads"]
    script:
        "../scripts/tophat_wrapper.py"

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

include: "bam_index.smk"
