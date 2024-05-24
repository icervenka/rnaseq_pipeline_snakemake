def get_count_output_files(wildcards):
    return [rules.count.output.counts, rules.counts_to_matrix.output]


def get_count_log_files(wildcards):
    return rules.move_count_summary.output


stranded_str = [
    featurecounts_stranded(
        Metadata.query("sample == @x").stranded.dropna().unique()[0]
    )
    for x in Samples
]

stranded_str = ",".join(stranded_str)

# TODO add stranded_str

rule count:
    input:
        bam=expand(rules.align_out.output, sample=Samples),
        gtf=config["gtf"],
    output:
        counts=opj(COUNT_OUTDIR + FEATURECOUNTS_COUNT_NAME),
        summary=opj(COUNT_OUTDIR + FEATURECOUNTS_SUMMARY_NAME),
    params:
        metadata=Metadata,
        # stranded="A"
        standard=featurecounts_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"]),
    log:
        opj(COUNT_LOG_OUTDIR + FEATURECOUNTS_LOG_FILES[0]),
    threads: 
        config["threads"]
    conda:
        CONDA_COUNT_GENERAL_ENV
    script:
        "../scripts/featurecounts_wrapper.py"


rule counts_to_matrix:
    input:
        rules.count.output.counts,
    output:
        opj(COUNT_OUTDIR + COMMON_COUNT_NAME)
    shell:
        # | awk '{{gsub("bam/","",$0); print}}'
        "cat {input} | head -n -1 | cut -f 1,7- > {output}"


rule move_count_summary:
    input:
        rules.count.output.summary,
    output:
        opj(COUNT_LOG_OUTDIR + FEATURECOUNTS_SUMMARY_NAME)
    shell:
        "mv {input} {output}"
