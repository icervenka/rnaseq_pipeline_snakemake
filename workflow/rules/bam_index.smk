def get_bam_index_files(wildcards):
    return expand(rules.bam_index.output, sample=Samples)

rule bam_index:
    input:
        rules.align_out.output
    output:
        ALIGN_OUTDIR + "{sample}/"+ COMMON_BAM_NAME + ".bam.bai"
    threads:
        config['threads']
    conda:
       CONDA_SHARED_ENV
    shell:
        "samtools index -@ {threads} {input} {output}"
