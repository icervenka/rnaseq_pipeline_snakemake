# TODO strandedness --rf --fr

rule count:
    input:
        bam=rules.align.output.bam,
        gtf=config['gtf']
    output:
        COUNT_OUTDIR + "{sample}/transcripts.gtf"
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
            "-o {output} "
            "{input.bam} "
        )

rule assemble:
    input:
        expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples)
    output:
        COUNT_OUTDIR + MERGE_OUTDIR + "gtf_assembly.txt"
    shell:
        "printf '%s\n' {input} >> {output}"

rule merge:
    input:
		merged=rules.assemble.output,
		gtf=config['gtf']
    output:
		COUNT_OUTDIR + MERGE_OUTDIR + "merged.gtf"
    threads:
        config['threads']
    run:
        shell(
            "stringtie "
			"--merge "
            "-p {threads} "
            "-G {input.gtf} "
            "-o {output} "
            "{input.merged} "
        )

rule count_merged:
    input:
		bam=rules.align.output.bam,
		merged=rules.merge.output
    output:
		DIFFEXP_OUTDIR + DIFFEXP_ANALYSIS + "/{sample}/final.gtf"
    threads:
        config['threads']
    run:
        shell(
            "stringtie "
			"-e "
			"-B "
            "-p {threads} "
            "-G {input.merged} "
            "-o {output} "
            "{input.bam} "
        )

rule all_count:
    input:
        expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples),
		expand(COUNT_OUTDIR + "{sample}/final.gtf", sample=Samples),
        COUNT_OUTDIR + MERGE_OUTDIR + "gtf_assembly.txt",
        COUNT_OUTDIR + MERGE_OUTDIR + "merged.gtf"
    output:
        touch(LOG_DIR + "count.completed")
