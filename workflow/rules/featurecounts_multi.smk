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
        counts=COUNT_OUTDIR + "gene_counts.txt",
        summary=COUNT_OUTDIR + "gene_counts.txt.summary",
    params:
        metadata=Metadata,
        # stranded="A"
        standard=featurecounts_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"]),
    log:
        COUNT_LOG_OUTDIR + "featurecounts_log.txt",
    threads: 
        config["threads"]
    conda:
        CONDA_COUNT_GENERAL_ENV
    script:
        "../scripts/featurecounts_wrapper.py"


rule counts_to_matrix:
    input:
        COUNT_OUTDIR + "gene_counts.txt",
    output:
        COUNT_OUTDIR + "counts.tsv",
    shell:
        # | awk '{{gsub("bam/","",$0); print}}'
        "cat {input} | head -n -1 | cut -f 1,7- > {output}"


rule move_count_summary:
    input:
        COUNT_OUTDIR + "gene_counts.txt.summary",
    output:
        COUNT_LOG_OUTDIR + "gene_counts.txt.summary",
    shell:
        "mv {input} {output}"
