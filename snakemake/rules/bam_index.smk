def get_bam_index_files(wildcards):
    return expand(ALIGN_OUTDIR + "{sample}/"+ COMMON_BAM_NAME + ".bam.bai", sample=Samples)

rule bam_index:
    input:
        rules.rename_bam.output
    output:
        ALIGN_OUTDIR + "{sample}/"+ COMMON_BAM_NAME + ".bam.bai"
    threads:
        config['threads']
    shell:
        "samtools index -@ {threads} {input} {output}"
