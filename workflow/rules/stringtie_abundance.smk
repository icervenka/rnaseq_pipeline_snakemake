

def get_count_output_files(wildcards):
    return (
        expand(rules.count.output.counts, sample=Samples) + 
        expand(rules.count.output.gtf, sample=Samples) + 
        expand(rules.count.output.ballgown, sample=Samples) + 
        get_stringtie_processed_output_files(wildcards)
    )

def get_count_log_files(wildcards):
    return []

rule count:
    input:
        bam=rules.align_out.output,
        gtf=config["gtf"]
    output:
        counts=opj(COUNT_OUTDIR, "{sample}", STRINGTIE_COUNT_NAME),
        gtf=opj(COUNT_OUTDIR, "{sample}", STRINGTIE_GTF_FILE),
        ballgown=expand(opj(COUNT_OUTDIR, "{{sample}}", "{file}"), file=BALLGOWN_INPUT_FILES)
    params:
        stranded=lambda wildcards: stranded_param(wildcards, "stringtie"),
        standard=stringtie_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    threads:
        config["threads"]
    conda:
        CONDA_COUNT_GENERAL_ENV
    shell:
        """
        stringtie \
        {params.stranded} \
        {params.standard} \
        {params.extra} \
        -p {threads} \
        -G {input.gtf} \
        -B \
        -e \
        -A {output.counts}\
        -o {output.gtf} \
        {input.bam}
        """


include: "stringtie_process_counts.smk"


