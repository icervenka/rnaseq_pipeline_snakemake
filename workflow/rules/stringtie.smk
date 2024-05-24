# TODO extra params to stringie can only go to one rule. If I want to pass more
## extra parameters it has to be reflected in the config
# TODO gather data and counts to matrix might go to separate rule since they are
## the same for each type of stringtie
# TODO stringtie guided and denovo only differ in -G in the first step, find
## a way to do it with param

def get_count_output_files(wildcards):
    return (
        get_stringtie_processed_output_files(wildcards) + 
        expand(rules.count.output.counts, sample=Samples) + 
        expand(rules.count.output.gtf, sample=Samples) + 
        expand(rules.count.output.ballgown, sample=Samples)
    )

def get_count_log_files(wildcards):
    return []

rule count:
    input:
        bam=rules.align_out.output,
        gtf=config["gtf"]
    output:
        counts=COUNT_OUTDIR + "{sample}/" + STRINGTIE_COUNT_NAME,
        gtf=COUNT_OUTDIR + "{sample}/" + STRINGTIE_GTF_FILE,
        ballgown=expand(COUNT_OUTDIR + "{{sample}}/{file}", file=BALLGOWN_INPUT_FILES)
    params:
        #stranded=stringtie_stranded,
        #standard=stringtie_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    threads:
        config["threads"]
    conda:
        CONDA_COUNT_GENERAL_ENV
    shell:
        """
        stringtie \
        {params.extra} \
        -p {threads} \
        -G {input.gtf} \
        -B \
        -e \
        -A {output.counts}\
        -o {output.gtf} \
        {input.bam} 
        """
            # "{params.stranded} "
			# "{params.standard} "


include "stringtie_process_counts.smk"


