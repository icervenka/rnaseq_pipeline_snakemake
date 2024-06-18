# requires that the directory with genomic fasta is writable
rule fasta_index:
    input:
        config['fasta']
    output:
        config['fasta'] + ".fai"
    conda:
        CONDA_SHARED_ENV
    shell:
        "samtools faidx -o {output} {input}"

# load subsampling rules if subsampling is specified
if is_set_subsample(config['preprocess']['subsample']):
    include: "subsample.smk"
    # change the current fastq dir for further processing
    FASTQ_CURRENT_DIR = FASTQ_PREPROCESSED_DIR
else:
    def get_subsample_pe_output_files(wildcards):
        return []

    def get_subsample_se_output_files(wildcards):
        return []

    def get_subsample_log_files(wildcards):
        return []

# load trimming rules if trimmer is specified
if is_set_trimmer(config['trim']['trimmer']):
    include: config['trim']['trimmer'] + ".smk"
    # change the current fastq dir for further processing
    FASTQ_CURRENT_DIR = FASTQ_TRIMMED_DIR
else:
    def get_trim_pe_output_files(wildcards):
        return []

    def get_trim_se_output_files(wildcards):
        return []

    def get_trim_log_files(wildcards):
        return []