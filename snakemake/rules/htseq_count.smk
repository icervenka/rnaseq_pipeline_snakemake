def get_count_output_files(wildcards):
    return [COUNT_LOG_OUTDIR + "htseq_count.log", COUNT_OUTDIR + "counts.tsv"]

def get_count_log_files(wildcards):
    return COUNT_LOG_OUTDIR + "gene_counts.txt.summary"

rule count:
    input:
        bam=rules.rename_bam.output,
        gtf=config["gtf"]
    output:
        temp(COUNT_OUTDIR + "{sample}_counts.txt")
    params:
        extra=htseq_params
    threads:
        config["threads"]
    run:
        stranded = Metadata.query('sample == @wildcards.sample')stranded.dropna().unique()[0]

        shell(
            "htseq-count "
            "-n {threads} "
            "-f bam "
            "-r pos "
            "-s {stranded} "
            "{params.extra} "
            "{input.bam} "
            "{input.gtf} "
            "> {output}
        )

rule counts_to_matrix:
    input:
        expand(rules.count.output, sample=Samples)
    output:
        COUNT_OUTDIR + "counts.tsv"
    log:
        COUNT_LOG_OUTDIR + "htseq_count.log"
    script:
        "../scripts/htseq_count_gather.py"
