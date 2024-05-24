def get_count_output_files(wildcards):
    return (
        [#rules.counts_to_matrix.output.counts_gene,
        #rules.counts_to_matrix.output.counts_transcript, 
        rules.gather_data.output.tpm,
        rules.gather_data.output.fpkm,
        rules.gather_data.output.samples_combined,
        rules.merge.output] +
        expand(rules.count.output.counts, sample=Samples) + 
        expand(rules.count.output.gtf, sample=Samples) + 
        expand(rules.count.output.ballgown, sample=Samples)
    )


def get_count_log_files(wildcards):
    return []


rule assemble_transcripts:
    input:
        bam=rules.align.output.bam,
        gtf=config['gtf']
    output:
        gtf=COUNT_OUTDIR + "{sample}/" + STRINGTIE_GTF_FILE,
	params:
        #stranded=stringtie_stranded,
        #standard=stringtie_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    threads:
        config['threads']
    conda:
        CONDA_COUNT_GENERAL_ENV
    shell:
        """
        stringtie \
        -p {threads} \
        -l {wildcards.sample} \
        -o {output.gtf} \
        {input.bam}
        """

        # {params.stranded} \
        # {params.extra} \"


rule merge:
    input:
		samples=expand(rules.assemble_transcripts.output, sample=Samples)
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
            "-o {output} "
            "{input.samples} "
        )


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
