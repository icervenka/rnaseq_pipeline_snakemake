# requires that the directory with genomic fasta is writable
rule fasta_index:
    input:
        config['fasta']
    output:
        config['fasta'] + ".fai"
    shell:
        "samtools faidx -o {output} {input}"
