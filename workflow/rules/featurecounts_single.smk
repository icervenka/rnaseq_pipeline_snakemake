def get_count_output_files(wildcards):
    return expand(rules.count.output.counts, sample=Samples) + [
        rules.counts_to_matrix.output.counts_tsv]


def get_count_log_files(wildcards):
    return (
        expand(rules.count.output.summary, sample=Samples) + 
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
        extra=featurecounts_params,
        # stranded="A",
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
        counts=expand(rules.count.output.counts, sample=Samples),
    output:
        counts_tsv=COUNT_OUTDIR + "counts.tsv",
    shell:
        "touch {output}"


# TODO think about merging into one file
rule move_count_summary:
    input:
          summary=rules.count.output.summary,
    output:
         COUNT_LOG_OUTDIR + "{sample}/" + "gene_counts.txt.summary"
    shell:
        "mv {input} {output}"
