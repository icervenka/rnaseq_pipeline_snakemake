def get_count_output_files(wildcards):
    return expand(rules.count.output, sample=Samples) + rules.counts_to_matrix.output


def get_count_log_files(wildcards):
    return expand(rules.move_count_log.output, sample=Samples)


rule count:
    input:
        bam=rules.align_out.output,
        gtf=config["gtf"]
    output:
        COUNT_OUTDIR + "{sample}/gene_counts.txt"
    params:
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
        {params.standard} \
        {params.extra} \
        {input.bam} \
        {input.gtf} \
        > {output}
        """


rule move_count_log:
    input:
        rules.count.output
    output:
        COUNT_LOG_OUTDIR + "{sample}/" + "htseq.log"
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
        COUNT_OUTDIR + "counts.tsv"
    params:
        kind="count"
    conda:
        CONDA_R_GENERAL_ENV
    script:
        "../scripts/htseq_gather.R"

