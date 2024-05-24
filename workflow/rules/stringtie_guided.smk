def get_count_output_files(wildcards):
    return (
        # [rules.counts_to_matrix.output.counts_gene,
        #rules.counts_to_matrix.output.counts_transcript, 
        # rules.gather_data.output.tpm,
        # rules.gather_data.output.fpkm,
        # rules.gather_data.output.samples_combined,
        get_stringtie_processed_output_files(wildcards) +
        [rules.merge.output] +
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
        -G {input.gtf} \
        -o {output.gtf} \
        {input.bam}
        """

        # {params.stranded} \
        # {params.extra} \"


rule merge:
    input:
		samples=expand(rules.assemble_transcripts.output, sample=Samples)
		gtf=config['gtf']
    output:
		COUNT_OUTDIR + STRINGTIE_MERGED_FILE
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
        --merge \
        -p {threads} \
        -G {input.gtf} \
        -o {output} \
        {input.samples}
        """

        # {params.stranded} \
        # {params.extra} \


rule count:
    input:
		bam=rules.align.output.bam,
		merged=rules.merge.output
    output:
        counts=COUNT_OUTDIR + "{sample}/" + STRINGTIE_COUNT_NAME,
		gtf=COUNT_OUTDIR + "{sample}/" + STRINGTIE_GTF_FILE
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
        -B \
        -e \
        -p {threads} \
        -G {input.merged} \
        -A {output.counts} \
        -o {output} \
        {input.bam}
        """

        # {params.stranded} \
        # {params.extra} \


include "stringtie_process_counts.smk"
