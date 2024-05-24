def get_count_output_files(wildcards):
    return (
        get_stringtie_processed_output_files(wildcards) +
        expand(rules.count.output.counts, sample=Samples) + 
        expand(rules.count.output.gtf, sample=Samples) + 
        expand(rules.count.output.ballgown, sample=Samples) +
        rules.merge.output
    )


def get_count_log_files(wildcards):
    return []


rule assemble_transcripts:
    input:
        bam=rules.align_out.output,
        gtf=config["gtf"]
    output:
        opj(COUNT_OUTDIR, "{sample}", "temp_" + STRINGTIE_GTF_FILE),
    params:
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
        -o {output} \
        {input.bam} 
        """


rule merge:
    input:
        sample_gtfs=expand(rules.assemble_transcripts.output, sample=Samples),
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
    shell: ":"
        """
        stringtie \
        --merge \
        -p {threads} \
        -G {input.gtf} \
        -o {output} \
        {input.sample_gtfs}
        """

        # {params.stranded} \
        # {params.extra} \


rule count:
    input:
        bam=rules.align_out.output,
        merged=rules.merge.output
    output:
        counts=opj(COUNT_OUTDIR, "{sample}", STRINGTIE_COUNT_NAME),
        gtf=opj(COUNT_OUTDIR, "{sample}", STRINGTIE_GTF_FILE),
        ballgown=expand(opj(COUNT_OUTDIR, "{{sample}}", "{file}"), file=BALLGOWN_INPUT_FILES)
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
        -o {output.gtf} \
        {input.bam}
        """

    # {params.stranded} \
    # {params.extra} \


include: "stringtie_process_counts.smk"

