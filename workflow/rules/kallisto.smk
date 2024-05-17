def get_align_output_files(wildcards):
    return (
        expand(
            ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".h5", sample=Samples)
        + expand(
            ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".tsv", sample=Samples)
        + expand(
            ALIGN_OUTDIR + "{sample}/" + "pseudoalignments.bam", sample=Samples)
    )


def get_align_log_files(wildcards):
    return expand(ALIGN_LOG_OUTDIR + "{sample}/" + "run_info.json", sample=Samples)


def get_bam_index_files(wildcards):
    return expand(
        ALIGN_OUTDIR + "{sample}/" + "pseudoalignments.bam.bai", sample=Samples
    )

# TODO fix passing stranded info
# move log to logs

rule kallisto_quant:
    input:
        sample=get_fq,
        index=config["index"],
        gtf=config["gtf"],
    output:
        h5=ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".h5",
        tsv=ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".tsv",
        bam=ALIGN_OUTDIR + "{sample}/" + "pseudoalignments.bam",
        log=ALIGN_OUTDIR + "{sample}/" + "run_info.json",
    params:
        metadata=Metadata,
        fastq_dir=FASTQ_INPUT_DIR,
        outdir=ALIGN_OUTDIR + "{sample}/",
        fragment_info=config_extra["align"]["kallisto_single_fragment_info"],
        #stranded=kallisto_params,
        extra=has_extra_config(config["align"]["extra"], config_extra["align"])
    threads: 
        config["threads"]
    conda:
        CONDA_ALIGN_GENERAL_ENV
    script:
        "../scripts/kallisto_wrapper.py"



rule move_align_log:
    input:
        rules.kallisto_quant.output.log,
    output:
        ALIGN_LOG_OUTDIR + "{sample}/" + "run_info.json",
    shell:
        "mv {input} {output}"

rule bam_index:
    input:
        rules.kallisto_quant.output.bam
    output:
        ALIGN_OUTDIR + "{sample}/"+ "pseudoalignments.bam.bai"
    threads:
        config['threads']
    conda:
        CONDA_SHARED_ENV
    shell:
        "samtools index -@ {threads} {input} {output}"

# include: "bam_index.smk"