def get_count_output_files(wildcards):
    return  (
        expand(rules.cufflinks.output, sample=Samples) +
        rules.cuffnorm.output
    )


def get_count_log_files(wildcards):
    return (
        expand(rules.cufflinks.log, sample=Samples) +
        rules.cuffnorm.log
    )


rule cufflinks:
    input:
        bam=rules.align_out.output,
    output:
        opj(COUNT_OUTDIR, "{sample}", CUFFLINKS_GTF_FILE)
    params:
        outdir=COUNT_OUTDIR,
        gtf=config["gtf"],
        fasta=config["fasta"],
        stranded=lambda wildcards: stranded_param(wildcards, "cufflinks"),
        standard=cufflinks_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    log:
        expand(opj(COUNT_LOG_OUTDIR, "{{sample}}", "{log}"), log=CUFFLINKS_LOG_FILES)
    threads:
        config["threads"]
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    shell:
        """
        cufflinks \
        -q \
        --no-update-check \
        --no-length-correction --no-effective-length-correction \
        -p {threads} \
        --compatible-hits-norm \
        --frag-bias-correct {params.fasta} \
        --library-type {params.stranded} \
        -G {params.gtf} \
        {params.extra} \
        -o {params.outdir}{wildcards.sample} \
        {input.bam} \
        > {log} 2>&1
        """


rule cuffnorm:
    input:
        bam=expand(rules.align_out.output, sample=Samples),
        gtf=config['gtf']
    output:
        expand(opj(COUNT_OUTDIR, "cuffnorm", "{files}"), files=CUFFNORM_COUNT_NAMES) 
    params:
        outdir=opj(COUNT_OUTDIR, "cuffnorm"),
        metadata=Metadata,
        labels=",".join(Samples),
        extra=config_extra["count_other_rules"]["cuffnorm_extra"]
    log:
        expand(opj(COUNT_LOG_OUTDIR, "{log}"), log=CUFFNORM_LOG_FILES)
    threads:
        config['threads']
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    script:
        # used a wrapper to unify strandedness info
        "../scripts/cuffnorm_wrapper.py"


# include: "counts_to_matrix.smk"