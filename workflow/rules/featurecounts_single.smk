def get_count_output_files(wildcards):
    return expand(rules.count.output.counts, sample=Samples) + [
        rules.counts_to_matrix.output]


def get_count_log_files(wildcards):
    return (
        expand(rules.move_count_summary.output, sample=Samples) + 
        expand(rules.count.log, sample=Samples)
    )

# TODO fix stranded

rule count:
    input:
        bam=rules.align_out.output,
        gtf=config["gtf"],
    output:
        counts=COUNT_OUTDIR + "{sample}/" + "gene_counts.txt",
        summary=COUNT_OUTDIR + "{sample}/" + "gene_counts.txt.summary",
    params:
        # stranded="A"
        standard=featurecounts_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    log:
        COUNT_LOG_OUTDIR + "{sample}/" + "featurecounts_log.txt", 
    threads: 
        config["threads"]
    conda:
        CONDA_COUNT_GENERAL_ENV
    script:
        "../scripts/featurecounts_wrapper.py"

# TODO 
rule counts_to_matrix:
    input:
        expand(rules.count.output.counts, sample=Samples),
    output:
        COUNT_OUTDIR + "counts.tsv",
    shell:
        "touch {output}"


# TODO think about merging into one file
rule move_count_summary:
    input:
        rules.count.output.summary,
    output:
        COUNT_LOG_OUTDIR + "{sample}/" + "gene_counts.txt.summary"
    shell:
        "mv {input} {output}"