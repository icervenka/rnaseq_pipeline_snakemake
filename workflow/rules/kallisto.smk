def get_align_output_files(wildcards):
    return  (
        expand(rules.align.output.h5,sample=Samples) + 
        expand(rules.align.output.tsv, sample=Samples) + 
        expand(rules.align_out.output, sample=Samples) +
        get_counts_to_matrix_output_files(wildcards)
    )


def get_align_log_files(wildcards):
    return expand(rules.move_align_log.output, sample=Samples)



rule align:
    input:
        sample=get_fq,
        gtf=config["gtf"],
    output:
        h5=opj(ALIGN_OUTDIR, "{sample}", KALLISTO_QUANT_FILE + ".h5"),
        tsv=opj(ALIGN_OUTDIR, "{sample}", KALLISTO_QUANT_FILE + ".tsv"),
        bam=opj(ALIGN_OUTDIR, "{sample}/", KALLISTO_BAM_FILE + ".bam"),
    params:
        metadata=Metadata,
        fastq_dir=FASTQ_CURRENT_DIR,
        index=config["index"],
        outdir=opj(ALIGN_OUTDIR, "{sample}"),
        fragment_info=config_extra["align"]["kallisto_single_fragment_info"],
        stranded=lambda wildcards: stranded_param(wildcards, "kallisto"),
        extra=has_extra_config(config["align"]["extra"], config_extra["align"])
    threads: 
        config["threads"]
    log:
        expand(opj(ALIGN_OUTDIR, "{{sample}}", "{logfiles}"), logfiles = KALLISTO_LOGFILES)
    conda:
        CONDA_ALIGN_GENERAL_ENV
    script:
        "../scripts/kallisto_wrapper.py"


rule move_align_log:
    input:
        rules.align.log
    output:
        expand(opj(ALIGN_LOG_OUTDIR, "{{sample}}", "{logfiles}"), logfiles = KALLISTO_LOGFILES)
    shell:
        "mv {input} {output}"


rule align_out:
    input:
        rules.align.output.bam
    output:
        opj(ALIGN_OUTDIR, "{sample}", COMMON_BAM_FILE + ".bam")
    shell:
        """
        mv {input} {output}
        """

include: "bam_index.smk"
include: "counts_to_matrix.smk"