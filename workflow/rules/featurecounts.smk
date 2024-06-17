def get_count_output_files(wildcards):
    return (
        expand(rules.count.output.counts, sample=Samples) + 
        rules.counts_to_matrix.output
    )


def get_count_log_files(wildcards):
    return (
        expand(rules.move_count_summary.output, sample=Samples) + 
        expand(rules.count.log, sample=Samples)
    )


rule count:
    input:
        bam=rules.align_out.output
    output:
        counts=opj(COUNT_OUTDIR, "{sample}", FEATURECOUNTS_COUNT_FILE),
        summary=opj(COUNT_OUTDIR, "{sample}", FEATURECOUNTS_SUMMARY_FILE)
    params:
        gtf=config["gtf"],
        fasta=config["fasta"],
        stranded=lambda wildcards: stranded_param(wildcards, "featurecounts"),
        standard=featurecounts_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    log:
        expand(opj(COUNT_LOG_OUTDIR, "{{sample}}", "{log}"), log=FEATURECOUNTS_LOG_FILES)
    threads: 
        config["threads"]
    conda:
        CONDA_COUNT_GENERAL_ENV
    shell:
        """
        featureCounts \
        -a {params.gtf} \
        -o {output.counts} \
        -T {threads} \
        -G {params.fasta} \
        -s {params.stranded} \
        {params.standard} \
        {params.extra} \
        {input.bam} \
        2> {snakemake.log}
        """


rule counts_to_matrix:
    input:
        expand(rules.count.output.counts, sample=Samples),
    output:
        opj(COUNT_OUTDIR, COMMON_COUNT_FILE),
    conda:
        CONDA_R_GENERAL_ENV
    script:
        "../scripts/featurecounts_gather.R"


rule move_count_summary:
    input:
        rules.count.output.summary,
    output:
        opj(COUNT_LOG_OUTDIR, "{sample}", FEATURECOUNTS_SUMMARY_FILE)
    shell:
        "mv {input} {output}"
