ruleorder: fastq_dump_pe > fastq_dump_se
# uses prefetch-validate combo because fasterq-dump has no way of
# checking if the file was downloaded correctly due to network errors

rule sra_download:
    input:
        ancient(config["metadata"])
    output:
        opj(SRA_DIR, "{sra}")
    params:
        get_sra
    threads:
        1
    conda:
        CONDA_SRA_TOOLS_ENV
    script:
        opj(CD2UP, WRAPPER_DIR, "sra_download_wrapper.py")
        

rule fastq_dump_se:
    input:
        get_single_end_sra
    output:
        opj(FASTQ_DIR, "{sra}.fastq")
    threads:
        config['threads']
    conda:
        CONDA_SRA_TOOLS_ENV
    shell:
        "fasterq-dump -e {threads} -s {input} -o {output}"

rule fastq_dump_pe:
    input:
        get_paired_end_sra
    output:
        expand(opj(FASTQ_DIR, "{{sra}}", "{read}.fastq"), read=['_1', "_2"])
    params:
        fastq_dir = FASTQ_DIR
    threads:
        config['threads']
    conda:
        CONDA_SRA_TOOLS_ENV
    shell:
        "fasterq-dump -e {threads} -3 {input} -O {params}"

rule compress_fastq:
    input:
        opj(FASTQ_DIR, "{file}.fastq")
    output:
        opj(FASTQ_DIR, "{file}.fastq.gz")
    threads:
        config['threads']
    conda:
        CONDA_SRA_TOOLS_ENV
    shell:
        "pigz --best -p {threads} {input}"