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
        bam=rules.align_out.output,
        gtf=config["gtf"],
    output:
        counts=opj(COUNT_OUTDIR, "{sample}", FEATURECOUNTS_COUNT_NAME),
        summary=opj(COUNT_OUTDIR, "{sample}", FEATURECOUNTS_SUMMARY_NAME)
    params:
        stranded=lambda wildcards: stranded_param(wildcards, "featurecounts"),
        standard=featurecounts_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    log:
        opj(COUNT_LOG_OUTDIR, "{sample}", FEATURECOUNTS_LOG_FILES[0]), 
    threads: 
        config["threads"]
    conda:
        CONDA_COUNT_GENERAL_ENV
    script:
        "../scripts/featurecounts_wrapper.py"

rule counts_to_matrix:
    input:
        expand(rules.count.output.counts, sample=Samples),
    output:
        opj(COUNT_OUTDIR, COMMON_COUNT_NAME),
    conda:
        CONDA_R_GENERAL_ENV
    script:
        "../scripts/featurecounts_gather.R"


# TODO think about merging into one file
rule move_count_summary:
    input:
        rules.count.output.summary,
    output:
        opj(COUNT_LOG_OUTDIR, "{sample}", FEATURECOUNTS_SUMMARY_NAME)
    shell:
        "mv {input} {output}"
