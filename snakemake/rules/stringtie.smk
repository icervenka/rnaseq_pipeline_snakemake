# TODO verify the ballgown filename
def get_count_output_files(wildcards):
    return expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples) +
        expand(COUNT_OUTDIR + "{sample}/gene_counts.tsv", sample=Samples) +
        expand(COUNT_OUTDIR + "{sample}/ballgown.ctab", sample=Samples)

def get_count_log_files(wildcards):
    return []

rule count:
    input:
        bam=rules.align.output.bam,
        gtf=config['gtf']
    output:
        gtf=COUNT_OUTDIR + "{sample}/transcripts.gtf",
        counts=COUNT_OUTDIR + "{sample}/gene_counts.tsv"
	params:
		extra=stringtie_params
    threads:
        config['threads']
    run:
        shell(
            "stringtie "
			"{params.extra} "
            "-p {threads} "
            "-G {input.gtf} "
			"-B "
            "-e "
            "-A {output.counts}"
            "-o {output.gtf} "
            "{input.bam} "
        )

rule counts_to_matrix:
    input:
        expand(rules.count.output.counts, sample=Samples)
    output:
        COUNT_OUTDIR + "counts.tsv"
    log:
        COUNT_LOG_OUTDIR + "stringtie_count.log"
    script:
        "../scripts/stringtie_count_gather.py"
