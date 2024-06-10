def get_count_output_files(wildcards):
    return (
        expand(rules.cufflinks.output, sample=Samples) +
        rules.cuffnorm.output +
        rules.cuffmerge.output
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
        denovo=cufflinks_denovo,
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
        {params.denovo} \
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
        denovo=cufflinks_denovo,
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
        {params.denovo} \
        {params.extra} \
        -o {params.outdir} \
        {input.gathered} \
        > {log.log} 2>&1
        """


# TODO hard coded file names
rule cuffmerge_move_log:
    input:
        rules.cuffmerge.log.runlog
    output:
        opj(COUNT_LOG_OUTDIR, "cuffmerge_run.log")
    params:
        outdir=opj(COUNT_OUTDIR, "logs")
    shell:
        """
        mv {input} {output}
        rm -rf {params.outdir}
        """

# rule cuffcompare:
#     input:
#         expand(rules.cufflinks.output, sample=Samples)
#     output:
#         COUNT_OUTDIR + "cuffcompare/" + "cuffcmp.loci"
#     params:
#         outdir=COUNT_OUTDIR + "cuffcompare/",
#         gtf=config["gtf"]
#     threads:
#         config['threads']
#     conda:
#         CONDA_COUNT_CUFFLINKS_ENV
#     shell:
#         "cuffcompare {input} -o {output}"


# NOTE cuffquant is only intermediate step and can be skipped
# rule cuffquant:
#     input:
#         bam=rules.align_out.output,
#         gtf=rules.cuffmerge.output.merged_gtf
#     output:
#         opj(COUNT_OUTDIR, "{sample}", CUFFQUANT_COUNT_NAME)
#     params:
#         outdir=COUNT_OUTDIR,
#         fasta=config["fasta"],
#         stranded=lambda wildcards: stranded_param(wildcards, "cufflinks"),
#         standard=cuffquant_params,
#         extra=config_extra["count_other_rules"]["cuffquant_extra"]
#     log:
#         expand(opj(COUNT_LOG_OUTDIR, "{{sample}}", "{log}"), log=CUFFQUANT_LOG_FILES)
#     threads:
#         config['threads']
#     conda:
#         CONDA_COUNT_CUFFLINKS_ENV
#     shell:
#         """
#         cuffquant \
#         --no-update-check \
#         -v \
#         -p {threads} \
#         --frag-bias-correct {params.fasta} \
#         --library-type {params.stranded} \
#         {params.standard} \
#         {params.extra} \
#         -o {params.outdir}{wildcards.sample}/ \
#         {input.gtf} \
#         {input.bam} \
#         > {log} 2>&1
#         """


# TODO only takes one stranded parameter, all files in the batch have to have same strandedness
rule cuffnorm:
    input:
        bam=expand(rules.align_out.output, sample=Samples),
        gtf=rules.cuffmerge.output.merged_gtf
    output:
        expand(opj(CUFFNORM_OUTDIR, "{files}"), files=CUFFNORM_COUNT_NAMES) 
    params:
        outdir=opj(CUFFNORM_OUTDIR),
        labels=",".join(Samples),
        # stranded=lambda wildcards: stranded_param(wildcards, "cufflinks"),
        extra=config_extra["count_other_rules"]["cuffnorm_extra"]
    log:
        expand(opj(COUNT_LOG_OUTDIR, "{log}"), log=CUFFNORM_LOG_FILES)
    threads:
        config['threads']
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    shell:
        """
        cuffnorm \
        -q \
        --no-update-check \
        --output-format simple-table \
        --library-norm-method classic-fpkm \
        -p {threads} \
        -o {params.outdir} \
        -L {params.labels} \
        {params.extra} \
        {input.gtf} \
        {input.bam} \
        > {log} 2>&1; 
        rm {params.outdir}/run.info
        """

# --library-type {params.stranded} \




include: "cuffnorm_to_raw.smk"

