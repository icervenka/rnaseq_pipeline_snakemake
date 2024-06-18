def get_count_output_files(wildcards):
    return (
        expand(rules.cufflinks.output, sample=Samples) +
        rules.cuffmerge.output + 
        rules.cuffnorm.output +
        get_cuffnorm_to_matrix(wildcards)
    )


def get_count_log_files(wildcards):
    return (
        expand(rules.cufflinks.log, sample=Samples) +
        rules.cuffnorm.log +
        rules.cuffmerge.log.log +
        rules.cuffmerge_move_log.output
    )


rule cufflinks:
    input:
        bam=rules.align_out.output
    output:
        opj(COUNT_OUTDIR, "{sample}", CUFFLINKS_GTF_FILE)
    params:
        outdir=COUNT_OUTDIR,
        fasta=config["fasta"],
        stranded=lambda wildcards: stranded_param(wildcards, "cufflinks"),
        # supplies gtf file
        assembly=cufflinks_assembly,
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
        -p {threads} \
        --library-type {params.stranded} \
        {params.assembly} \
        {params.standard} \
        {params.extra} \
        -o {params.outdir}{wildcards.sample} \
        {input.bam} \
        > {log} 2>&1
        """

#         --frag-bias-correct {params.fasta} \

rule gather_gtf:
    input:
        expand(rules.cufflinks.output, sample=Samples)
    output:
        temp(COUNT_LOG_OUTDIR + "gtf_assembly.txt")
    shell:
        "printf '%s\n' {input} >> {output}"


rule cuffmerge:
    input:
        gathered=rules.gather_gtf.output,
    output:
        merged_gtf=opj(COUNT_OUTDIR, CUFFLINKS_MERGED_FILE)
    params:
        outdir=COUNT_OUTDIR,
        fasta=config["fasta"],
        gtf=config["gtf"],
        # assembly=cufflinks_assembly,
        extra=config_extra["count_other_rules"]["cuffmerge_extra"]
    log:
        runlog=opj(COUNT_OUTDIR, "logs", "run.log"),
        log=expand(opj(COUNT_LOG_OUTDIR, "{log}"), log=CUFFMERGE_LOG_FILES)
    threads:
        config['threads']
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    shell:
        """
        cuffmerge \
        -p {threads} \
        -s {params.fasta} \
        -g {params.gtf} \
        {params.extra} \
        -o {params.outdir} \
        {input.gathered} \
        > {log.log} 2>&1
        """


rule cuffmerge_move_log:
    input:
        rules.cuffmerge.log.runlog
    output:
        expand(opj(COUNT_LOG_OUTDIR, "{log}"), log=CUFFMERGE_RUN_LOG_FILES)
    params:
        outdir=opj(COUNT_OUTDIR, "logs")
    shell:
        """
        mv {input} {output}
        rm -rf {params.outdir}
        """

rule cuffcompare:
    input:
        rules.cuffmerge.output.merged_gtf
    output:
        gtf=opj(CUFFCOMPARE_OUTDIR, "cuffcmp" + CUFFCOMPARE_GTF_FILE),
        other=expand(opj(CUFFCOMPARE_OUTDIR, "cuffcmp" + "{files}"), files=CUFFCOMPARE_FILES)
    params:
        outprefix=opj(CUFFCOMPARE_OUTDIR, "cuffcmp"),
        gtf=config["gtf"]
    threads:
        config['threads']
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    shell:
        "cuffcompare {input} {params.gtf} -o {outprefix}"


rule cuffnorm:
    input:
        bam=expand(rules.align_out.output, sample=Samples),
        gtf=rules.cuffmerge.output.merged_gtf
    output:
        expand(opj(CUFFNORM_OUTDIR, "{files}"), files=CUFFNORM_COUNT_FILES) 
    params:
        outdir=opj(CUFFNORM_OUTDIR),
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
        opj(CD2UP, SCRIPT_DIR, "cuffnorm_wrapper.py")


include: "cuffnorm_to_matrix.smk"

