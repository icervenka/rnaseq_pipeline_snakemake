def get_count_output_files(wildcards):
    return [COUNT_OUTDIR + "gene_counts.txt", COUNT_OUTDIR + "counts.tsv"]

def get_count_log_files(wildcards):
    return COUNT_LOG_OUTDIR + "gene_counts.txt.summary"

rule count:
    input:
        bam=expand(rules.rename_bam.output, sample=Samples),
        gtf=config["gtf"]
    output:
        counts=COUNT_OUTDIR + "gene_counts.txt",
        summary=COUNT_OUTDIR + "gene_counts.txt.summary"
    params:
        extra=featurecounts_params
    log:
        COUNT_LOG_OUTDIR + "featurecounts_log.txt"
    threads:
        config["threads"]
    run:
        stranded_str = [ featurecounts_stranded(Metadata.query('sample == @x').
            stranded.dropna().unique()[0]) for x in Samples ]
        stranded_str = ','.join(stranded_str)
        shell(
            "featureCounts "
            "-a {input.gtf} "
            "-o {output.counts} "
            "-T {threads} "
            "-s {stranded_str} "
            "{params.extra} "
            "{input.bam} "
            "2> {log}"
        )

rule counts_to_matrix:
    input:
        counts=COUNT_OUTDIR + "gene_counts.txt"
    output:
        counts_tsv=COUNT_OUTDIR + "counts.tsv"
    shell:
        # | awk '{{gsub("bam/","",$0); print}}'
        "cat {input.counts} | head -n -1 | cut -f 1,7- > {output.counts_tsv}"

rule move_count_log:
    input:
        COUNT_OUTDIR + "gene_counts.txt.summary"
    output:
        COUNT_LOG_OUTDIR + "gene_counts.txt.summary"
    shell:
        "mv {input} {output}"
