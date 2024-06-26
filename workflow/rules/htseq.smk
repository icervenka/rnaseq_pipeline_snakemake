def get_count_output_files(wildcards):
    return expand(rules.count.output, sample=Samples) + rules.counts_to_matrix.output


def get_count_log_files(wildcards):
    return expand(rules.move_count_log.output, sample=Samples)


# TODO redirect the other out to log file
rule count:
    input:
        bam=rules.align_out.output,
    output:
        opj(COUNT_OUTDIR, "{sample}", HTSEQ_COUNT_FILE)
    params:
        gtf=config["gtf"],
        stranded=lambda wildcards: stranded_param(wildcards, "htseq"),
        standard=htseq_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    threads:
        config["threads"]
    conda:
        CONDA_COUNT_GENERAL_ENV
    shell:
        """
        htseq-count \
        -n {threads} \
        -f bam \
        -r pos \
        -s {params.stranded} \
        {params.standard} \
        {params.extra} \
        {input.bam} \
        {params.gtf} \
        > {output}
        """


rule move_count_log:
    input:
        rules.count.output
    output:
        expand(opj(COUNT_LOG_OUTDIR, "{{sample}}", "{log}"), log=HTSEQ_LOG_FILES)
    params:
        kind="log"
    conda:
        CONDA_R_GENERAL_ENV
    script:
        "../scripts/htseq_gather.R"


rule counts_to_matrix:
    input:
        expand(rules.count.output, sample=Samples)
    output:
        opj(COUNT_OUTDIR, COMMON_COUNT_FILE)
    params:
        kind="count"
    conda:
        CONDA_R_GENERAL_ENV
    script:
        "../scripts/htseq_gather.R"

