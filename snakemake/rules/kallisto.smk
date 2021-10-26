def get_align_output_files(wildcards):
    return
        expand(ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".h5", sample = Samples) +
        expand(ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".tsv", sample = Samples) +
        expand(ALIGN_OUTDIR + "{sample}/" + "pseudoalignments.bam", sample = Samples) +
        expand(ALIGN_OUTDIR + "{sample}/" + "pseudoalignments.bam.bai", sample = Samples)

def get_align_log_files(wildcards):
    return expand(ALIGN_LOG_OUTDIR + "{sample}/" + "run_info.json", sample = Samples)

rule kallisto_quant:
    input:
        sample=get_fq,
        index=config['index'],
        gtf=config['gtf']
    output:
        h5=ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".h5",
        tsv=ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".tsv",
        log=ALIGN_OUTDIR + "{sample}/" + "run_info.json"
    params:
        single_extra=config_extra['align']['kallisto_single_extra'],
        extra=kallisto_params
        # stranded=kallisto_stranded,
        # length_mean=config_extra['align']['kallisto_extra']['read_length_mean'],
        # length_variance=config_extra['align']['kallisto_extra']['read_length_variance'],
        # extra=config_extra['align']['kallisto_extra']['extra']
    threads:
        config['threads']
    script:
        "../scripts/kallisto_wrapper.py"


rule move_align_log:
    input:
        rules.align.output.log
    output:
        ALIGN_LOG_OUTDIR + "{sample}/" + "run_info.json"
    shell:
        "cp {input} {output}"

# rule all_align:
#     input:
#         expand(ALIGN_OUTDIR+"{sample}/"+KALLISTO_BAM_NAME+".h5", sample = Samples),
#         expand(ALIGN_OUTDIR+"{sample}/"+KALLISTO_BAM_NAME+".tsv", sample = Samples)
#     output:
#         touch(LOG_DIR + "align.completed"),
#         touch(LOG_DIR + "count.completed")
