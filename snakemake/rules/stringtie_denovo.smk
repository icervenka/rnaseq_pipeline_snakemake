def get_count_output_files(wildcards):
    return expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples) +
        expand(COUNT_OUTDIR + "{sample}/gene_counts.gtf", sample=Samples) +
        [COUNT_OUTDIR + "merged.gtf"]

def get_count_log_files(wildcards):
    return []

rule assemble_transcripts:
    input:
        bam=rules.align.output.bam,
    output:
        COUNT_OUTDIR + "{sample}/transcripts.gtf"
	params:
		extra=stringtie_assemble_params
    threads:
        config['threads']
    run:
        shell(
            "stringtie "
			"{params.extra} "
            "-p {threads} "
            "-l {wildcards.sample} "
            "-o {output.gtf} "
            "{input.bam} "
        )

# TODO merge requires strandedness
rule merge:
    input:
		samples=expand(rules.assemble_transcripts.output, sample=Samples)
    output:
		COUNT_OUTDIR + "merged.gtf"
    params:
        extra=stringtie_merge_params
    threads:
        config['threads']
    run:
        shell(
            "stringtie "
			"--merge "
            "{params.extra} "
            "-p {threads} "
            "-o {output} "
            "{input.samples} "
        )

# TODO count requires strandedness
rule count:
    input:
		bam=rules.align.output.bam,
		merged=rules.merge.output
    output:
		gtf=COUNT_OUTDIR + "{sample}/transcripts.gtf",
        counts=COUNT_OUTDIR + "{sample}/gene_counts.tsv"
    params:
        extra=stringtie_count_params
    threads:
        config['threads']
    run:
        shell(
            "stringtie "
			"-B "
            "-e "
            "{params.extra} "
            "-p {threads} "
            "-G {input.merged} "
            "-A {output.counts}"
            "-o {output} "
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
