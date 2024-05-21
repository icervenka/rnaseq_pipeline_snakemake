def get_align_output_files(wildcards):
    return expand(ALIGN_OUTDIR + "{sample}/" + SALMON_QUANT_NAME + ".sf", sample=Samples)

def get_align_log_files(wildcards):
    return expand(ALIGN_LOG_OUTDIR + "{sample}/{log}",
        sample=Samples, log=SALMON_LOG_FILES+["salmon_quant.log"])

def get_bam_index_files(wildcards):
    return []

# TODO fix stranded
# TODO log is being created anyway

rule align:
    input:
        fq=get_fq,
        gtf=config["gtf"]
    output: 
        bam=ALIGN_OUTDIR + "{sample}/" + SALMON_QUANT_NAME + ".sf",
        runlog=expand(ALIGN_OUTDIR + "{{sample}}/" + "{runlog}", runlog=SALMON_LOG_FILES)
    params:
        index=config["index"],
        metadata=Metadata,
        fastq_dir=FASTQ_INPUT_DIR,
        outdir=ALIGN_OUTDIR + "{sample}/",
        fragment_info=config_extra["align"]["salmon_single_fragment_info"],
        extra=has_extra_config(config["align"]["extra"], config_extra["align"]),
        stranded="A"
    log:
        ALIGN_LOG_OUTDIR + "{sample}/salmon_quant.log"
    threads:
        config["threads"]
    conda:
        CONDA_ALIGN_GENERAL_ENV
    script:
        "../scripts/salmon_wrapper.py"

# rule align_out:
#     input:
#         rules.align.output.bam
#     output:
#         ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".bam"
#     shell:
#         "mv {input} {output}"


rule move_align_log:
    input:
        rules.align.output.runlog
    output:
        expand(ALIGN_LOG_OUTDIR + "{{sample}}/{log}", log=SALMON_LOG_FILES)
    params:
        outdir=ALIGN_LOG_OUTDIR + "{sample}/"
    shell:
        "mv {input} {params.outdir}"

