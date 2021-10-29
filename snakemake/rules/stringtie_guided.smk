def get_count_output_files(wildcards):
    return expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples) +
        expand(COUNT_OUTDIR + "{sample}/gene_counts.tsv", sample=Samples) +
        expand(COUNT_OUTDIR + "{sample}/{ballgown}.ctab", sample=Samples,
            ballgown=BALLGOWN_INPUT_FILES) +
        [COUNT_OUTDIR + "merged.gtf",
            COUNT_OUTDIR + "counts.tsv",
            COUNT_OUTDIR + "fpkm.tsv",
            COUNT_OUTDIR + "samples_combined.tsv"]

def get_count_log_files(wildcards):
    return []

rule assemble_transcripts:
    input:
        bam=rules.align.output.bam,
        gtf=config['gtf']
    output:
        gtf=COUNT_OUTDIR + "{sample}/transcripts.gtf"
	params:
		stranded=stringtie_stranded,
        extra=config_extra['count']['stringtie_assemble_extra']
    threads:
        config['threads']
    run:
        shell(
            "stringtie "
            "{params.stranded} "
			"{params.extra} "
            "-p {threads} "
            "-l {wildcards.sample} "
            "-G {input.gtf} "
            "-o {output.gtf} "
            "{input.bam} "
        )

# TODO merge requires strandedness
rule merge:
    input:
		samples=expand(rules.assemble_transcripts.output, sample=Samples)
		gtf=config['gtf']
    output:
		COUNT_OUTDIR + "merged.gtf"
    params:
		stranded=stringtie_stranded,
        extra=config_extra['count']['stringtie_merge_extra']
    threads:
        config['threads']
    run:
        shell(
            "stringtie "
			"--merge "
            "{params.stranded} "
            "{params.extra} "
            "-p {threads} "
            "-G {input.gtf} "
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
		stranded=stringtie_stranded,
        extra=config_extra['count']['stringtie_count_extra']
    threads:
        config['threads']
    run:
        shell(
            "stringtie "
			"-B "
            "-e "
            "{params.stranded} "
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
        counts=COUNT_OUTDIR + "counts.tsv",
        fpkm=COUNT_OUTDIR + "fpkm.tsv",
        samples_combined=COUNT_OUTDIR + "samples_combined.tsv"
    script:
        "../scripts/stringtie_count_gather.R"
