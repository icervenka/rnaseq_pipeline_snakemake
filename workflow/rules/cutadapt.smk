ruleorder: cutadapt_pe > cutadapt_se

def get_trim_pe_output_files(wildcards):
    return expand(rules.cutadapt_pe.output, 
        filename=list(get_metadata("filename_sans_read", 1)),
        ext=list(get_metadata("filename_full_ext", 1))
    )

def get_trim_se_output_files(wildcards):
    return expand(rules.cutadapt_se.output, 
        filename=list(get_metadata("filename_sans_read", 0)),
        ext=list(get_metadata("filename_full_ext", 0))
    )


def get_trim_log_files(wildcards):
    return expand(LOG_DIR + "trim/{filename}.txt", 
        filename=list(
            Metadata["filename_sans_read"].unique()
        )
    )


rule cutadapt_se:
    input:
        lambda wildcards: get_single_fq(wildcards, FASTQ_CURRENT_DIR)
    output:
        opj(FASTQ_TRIMMED_DIR, "{filename}{ext}")
    params:
        adapters=config['trim']['adapters_single'],
        quality=config['trim']['quality'],
        extra=config['trim']['extra']
    threads: 
        4
    log:
        opj(TRIM_LOG_OUTDIR, "{filename}{ext}.txt")
    conda:
        CONDA_SHARED_ENV
    shell:
        """
        cutadapt \
        {params.adapters} \
        {params.quality} \
        {params.extra} \
        -j {threads} \
        -o {output} \
        {input} \
        > {log} "
        """

rule cutadapt_pe:
    input:
        lambda wildcards: get_paired_fq(wildcards, FASTQ_CURRENT_DIR)
    output:
        expand(opj(FASTQ_TRIMMED_DIR, "{{filename}}{read}{{ext}}"), 
            read=config["paired_read_strings"]
        )
    params:
        adapters=config['trim']['adapters_paired'],
        quality=config['trim']['quality'],
        extra=config['trim']['extra']
    log:
        opj(TRIM_LOG_OUTDIR, "{filename}{ext}.txt")
    threads: 
        4
    conda:
        CONDA_SHARED_ENV
    script:
        "../scripts/cutadapt_pe_wrapper.py"

