def get_align_output_files(wildcards):
    return expand(rules.align.output.bam, sample=Samples)

def get_align_log_files(wildcards):
    return (
        expand(rules.move_align_log.output.quantlog, sample=Samples) +
        expand(rules.move_align_log.output.runlog, sample=Samples)
    )

def get_bam_index_files(wildcards):
    return []


rule align:
    input:
        fq=get_fq,
        gtf=config["gtf"]
    output: 
        bam=opj(ALIGN_OUTDIR, "{sample}", SALMON_QUANT_FILE)
    params:
        index=config["index"],
        metadata=Metadata,
        fastq_dir=FASTQ_CURRENT_DIR,
        outdir=opj(ALIGN_OUTDIR, "{sample}"),
        fragment_info=config_extra["align"]["salmon_single_fragment_info"],
        stranded=salmon_stranded,
        extra=has_extra_config(config["align"]["extra"], config_extra["align"]),
    log:
        quantlog=expand(opj(ALIGN_OUTDIR, "{{sample}}", "logs", "{log}"), log=["salmon_quant.log"]),
        runlog=expand(opj(ALIGN_OUTDIR, "{{sample}}", "{log}"), log=SALMON_LOG_FILES)
    threads:
        config["threads"]
    conda:
        CONDA_ALIGN_GENERAL_ENV
    script:
        opj(CD2UP, WRAPPER_DIR, "salmon_wrapper.py")


# TODO snakemake needs to modify files to avoid cyclic dependencies, it can't do noop
# later when feeding this to diffexp it might become a problem, renaming to common files would work
# rule align_out:
#     input:
#         rules.align.output.bam
#     output:
#         ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_FILE + ".bam"
#     shell:
#         "mv {input} {output}"

rule move_align_log:
    input:
        quantlog=rules.align.log.quantlog,
        runlog=rules.align.log.runlog
    output:
        quantlog=expand(opj(ALIGN_LOG_OUTDIR, "{{sample}}", "{log}"), log=["salmon_quant.log"]),
        runlog=expand(opj(ALIGN_LOG_OUTDIR, "{{sample}}", "{log}"), log=SALMON_LOG_FILES)
    params:
        logdir=opj(ALIGN_OUTDIR, "{sample}", "logs"),
        outdir=opj(ALIGN_LOG_OUTDIR, "{sample}")
    shell:
        """
        mv {input.runlog} {params.outdir}
        cp {input.quantlog} {params.outdir}
        rm -rf {params.logdir}
        """


