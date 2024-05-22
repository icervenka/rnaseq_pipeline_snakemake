def get_align_output_files(wildcards):
    return  (
        expand(rules.align.output.h5,sample=Samples) + 
        expand(rules.align.output.tsv, sample=Samples) + 
        expand(rules.align_out.output.bam, sample=Samples)
    )


def get_align_log_files(wildcards):
    return expand(rules.move_align_log.output, sample=Samples)


# TODO fix passing stranded info

rule align:
    input:
        sample=get_fq,
        index=config["index"],
        gtf=config["gtf"],
    output:
        h5=ALIGN_OUTDIR + "{sample}/" + KALLISTO_QUANT_NAME + ".h5",
        tsv=ALIGN_OUTDIR + "{sample}/" + KALLISTO_QUANT_NAME + ".tsv",
        bam=ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".bam",
    params:
        metadata=Metadata,
        fastq_dir=FASTQ_INPUT_DIR,
        outdir=ALIGN_OUTDIR + "{sample}/",
        fragment_info=config_extra["align"]["kallisto_single_fragment_info"],
        #stranded=kallisto_params,
        extra=has_extra_config(config["align"]["extra"], config_extra["align"])
    threads: 
        config["threads"]
    log:
        expand(ALIGN_OUTDIR + "{{sample}}/{logfiles}", logfiles = KALLISTO_LOGFILES)
    conda:
        CONDA_ALIGN_GENERAL_ENV
    script:
        "../scripts/kallisto_wrapper.py"

rule move_align_log:
    input:
        rules.align.log
    output:
        expand(ALIGN_LOG_OUTDIR + "{{sample}}/{logfiles}", logfiles = KALLISTO_LOGFILES)
    shell:
        "mv {input} {output}"


# TODO snakemake needs to modify files to avoid cyclic dependencies, it can't do noop
# later when feeding this to diffexp it might become a problem, renaming to common files would work
rule align_out:
    input:
        bam=rules.align.output.bam
    output:
        bam=ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".bam"
    shell:
        """
        mv {input.bam} {output.bam}
        """

include: "bam_index.smk"
