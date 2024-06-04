def get_count_output_files(wildcards):
    return (
        get_stringtie_processed_output_files(wildcards) +
        expand(rules.count.output.counts, sample=Samples) + 
        expand(rules.count.output.gtf, sample=Samples) + 
        expand(rules.count.output.ballgown, sample=Samples) +
        rules.merge.output
    )


def get_count_log_files(wildcards):
    return (
        expand(rules.assemble_transcripts.log, sample=Samples)
        expand(rules.merge.log, sample=Samples)
        expand(rules.count.log, sample=Samples)
    )


rule assemble_transcripts:
    input:
        bam=rules.align_out.output
    output:
        gtf=temp(opj(COUNT_OUTDIR, "{sample}", "temp_" + STRINGTIE_GTF_FILE))
    params:
        stranded=lambda wildcards: stranded_param(wildcards, "stringtie"),
        denovo=stringtie_denovo,
        standard=stringtie_params,
        extra=config_extra["count_other_rules"]["stringtie_assemble_extra"]
    log:
        expand(opj(COUNT_LOG_OUTDIR, "{{sample}}", "{log}"), log=STRINGTIE_LOG_FILES)
    threads:
        config['threads']
    conda:
       CONDA_COUNT_GENERAL_ENV
    shell:
        """
        stringtie \
        -v \
        {params.stranded} \
        {params.denovo} \
        {params.standard} \
        {params.extra} \
        -p {threads} \
        -l {wildcards.sample} \
        -o {output} \
        {input.bam} 
        """

rule merge:
    input:
        sample_gtfs=expand(rules.assemble_transcripts.output.gtf, sample=Samples),
    output:
        merged_gtf=opj(COUNT_OUTDIR, STRINGTIE_MERGED_FILE)
    params:
        denovo=stringtie_denovo,
        extra=config_extra["count_other_rules"]["stringtie_merge_extra"]
    log:
        expand(opj(COUNT_LOG_OUTDIR, "{{sample}}", "{log}"), log=STRINGTIE_MERGE_LOG_FILES)
    threads:
        config['threads']
    conda:
        CONDA_COUNT_GENERAL_ENV
    shell: 
        """
        stringtie \
        --merge \
        -v \
        {params.denovo} \
        {params.extra} \
        -p {threads} \
        -o {output} \
        {input.sample_gtfs}
        """

# to compare  transcripts with reference
# gffcompare -r $RNA_REF_GTF -o gffcompare stringtie_merged.gtf
# cat gffcompare.stats


rule count:
    input:
        bam=rules.align_out.output,
        merged=rules.merge.output.merged_gtf
    output:
        counts=opj(COUNT_OUTDIR, "{sample}", STRINGTIE_COUNT_NAME),
        gtf=opj(COUNT_OUTDIR, "{sample}", STRINGTIE_GTF_FILE),
        ballgown=expand(opj(COUNT_OUTDIR, "{{sample}}", "{file}"), file=BALLGOWN_INPUT_FILES)
    params:
        stranded=lambda wildcards: stranded_param(wildcards, "stringtie"),
        standard=stringtie_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    log:
        expand(opj(COUNT_LOG_OUTDIR, "{{sample}}", "{log}"), log=STRINGTIE_COUNT_LOG_FILES)
    threads:
        config['threads']   
    conda:
        CONDA_COUNT_GENERAL_ENV     
    shell:
        """
        stringtie \
        -v \
        {params.stranded} \
        {params.standard} \
        {params.extra} \
        -B \
        -e \
        -p {threads} \
        -G {input.merged} \
        -A {output.counts} \
        -o {output.gtf} \
        {input.bam}
        """


include: "stringtie_process_counts.smk"

