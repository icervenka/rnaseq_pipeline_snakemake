def get_bam_index_files(wildcards):
    expand(ALIGN_OUTDIR+"{sample}/"+COMMON_BAM_NAME+".bam.bai", sample=Samples)

rule bam_index:
    input:
        rules.rename_bam.output[0]
    output:
        ALIGN_OUTDIR + "{sample}/"+ COMMON_BAM_NAME + "bam.bai"
    threads:
        config['threads']
    shell:
        "samtools index -@ {threads} {input} {output}"
