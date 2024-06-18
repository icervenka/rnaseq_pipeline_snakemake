ruleorder: subsample_pe > subsample_se

def get_subsample_pe_output_files(wildcards):
    return expand(rules.subsample_pe.output, 
        filename=list(get_metadata("filename_sans_read", 1)),
        ext=list(get_metadata("filename_full_ext", 1))
    )


def get_subsample_se_output_files(wildcards):
    return expand(rules.subsample_se.output, 
        filename=list(get_metadata("filename_sans_read", 0)),
        ext=list(get_metadata("filename_full_ext", 0))
    )


# TODO fix subsample output log files
def get_subsample_log_files(wildcards):
    return []


rule subsample_se:
    input:
        lambda wildcards: get_single_fq(wildcards, FASTQ_CURRENT_DIR)
    output:
        opj(FASTQ_PREPROCESSED_DIR, "{filename}{ext}")
    params:
        proportion=get_subsample_proportion(config['preprocess']['subsample']),
        extra=config['preprocess']['extra']
    log:
        opj(SUBSAMPLE_LOG_OUTDIR, "{filename}{ext}.log")
    threads: 
        1
    conda:
        CONDA_SHARED_ENV
    shell:
       "fq subsample --r1-dst {output} {params.proportion} {input} > {log} 2>&1"

rule subsample_pe:
    input:
        lambda wildcards: get_paired_fq(wildcards, FASTQ_CURRENT_DIR)
    output:
        expand(opj(FASTQ_PREPROCESSED_DIR, "{{filename}}{read}{{ext}}"), 
            read=config["paired_read_strings"]
        )
    params:
        proportion=get_subsample_proportion(config['preprocess']['subsample']),
        extra=config['preprocess']['extra']
    log:
        opj(SUBSAMPLE_LOG_OUTDIR, "{filename}{ext}.log")
    threads: 
        1
    conda:
        CONDA_SHARED_ENV
    script:
        opj(CD2UP, WRAPPER_DIR, "subsample_wrapper.py")

