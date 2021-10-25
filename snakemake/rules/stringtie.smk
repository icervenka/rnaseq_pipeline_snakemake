def get_count_output_files(wildcards):
    return [expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples) +
            expand(COUNT_OUTDIR + "{sample}/final.gtf", sample=Samples) +
            [COUNT_OUTDIR + "merged.gtf", DIFFEXP_OUTDIR + DIFFEXP_ANALYSIS + "/{sample}/final.gtf"]]

def get_count_log_files(wildcards):
    return COUNT_OUTDIR + "gtf_assembly.txt",

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
        COUNT_OUTDIR + "gtf_assembly.txt"
    shell:
        "printf '%s\n' {input} >> {output}"

rule merge:
    input:
		merged=rules.assemble.output,
		gtf=config['gtf']
    output:
		COUNT_OUTDIR + "merged.gtf"
    params:
        extra=config_extra['count']['stringtie_merge_extra']
    threads:
        config['threads']
    run:
        shell(
            "stringtie "
			"--merge "
            "{params.extra} "
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
    params:
        extra=config_extra['count']['stringtie_merge_count_extra']
    threads:
        config['threads']
    run:
        shell(
            "stringtie "
			"-B "
            "{params.extra} "
            "-p {threads} "
            "-G {input.merged} "
            "-o {output} "
            "{input.bam} "
        )
