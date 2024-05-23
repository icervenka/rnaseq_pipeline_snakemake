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
else:
    pass

# load trimmin rules if trimmer is specified
if is_set_trimmer(config['trim']['trimmer']):
    include: config['trim']['trimmer'] + ".smk"
else:
    pass