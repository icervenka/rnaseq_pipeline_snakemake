def get_count_output_files(wildcards):
    return (
        expand(rules.count.output.counts, sample=Samples) + 
        expand(rules.count.output.gtf, sample=Samples) + 
        expand(rules.count.output.ballgown, sample=Samples) +
        rules.merge.output +
        get_stringtie_processed_output_files(wildcards) +
        get_counts_to_matrix_output_files(wildcards)
    )


def get_count_log_files(wildcards):
    return (
        expand(rules.assemble_transcripts.log, sample=Samples) +
        expand(rules.count.log, sample=Samples) +
        rules.merge.log
    )


rule assemble_transcripts:
    input:
        bam=rules.align_out.output
    output:
        gtf=temp(opj(COUNT_OUTDIR, "{sample}", "temp_" + STRINGTIE_GTF_FILE))
    params:
        stranded=lambda wildcards: stranded_param(wildcards, "stringtie"),
        assembly=stringtie_assembly,
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
        {params.stranded} \
        {params.assembly} \
        {params.standard} \
        {params.extra} \
        -p {threads} \
        -l {wildcards.sample} \
        -o {output.gtf} \
        {input.bam} \
        > {log} 2>&1
        """

rule merge:
    input:
        sample_gtfs=expand(rules.assemble_transcripts.output.gtf, sample=Samples),
    output:
        merged_gtf=opj(STRINGTIE_OUTDIR, STRINGTIE_MERGED_FILE)
    params:
        assembly=stringtie_assembly,
        extra=config_extra["count_other_rules"]["stringtie_merge_extra"]
    log:
        expand(opj(COUNT_LOG_OUTDIR, "{log}"), log=STRINGTIE_MERGE_LOG_FILES)
    threads:
        config['threads']
    conda:
        CONDA_COUNT_GENERAL_ENV
    shell: 
        """
        stringtie \
        --merge \
        {params.assembly} \
        {params.extra} \
        -p {threads} \
        -o {output} \
        {input.sample_gtfs}
        > {log} 2>&1
        """

# to compare  transcripts with reference
# gffcompare -r $RNA_REF_GTF -o gffcompare stringtie_merged.gtf
# cat gffcompare.stats


rule count:
    input:
        bam=rules.align_out.output,
        merged=rules.merge.output.merged_gtf
    output:
        counts=opj(COUNT_OUTDIR, "{sample}", STRINGTIE_COUNT_FILE),
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
        {input.bam} \
        > {log} 2>&1
        """

include: "counts_to_matrix.smk"
include: "stringtie_process_counts.smk"