rule bam_index:
    input:
        rules.rename_bam.output[0]
    output:
        ALIGN_OUTDIR + "{sample}/"+ COMMON_BAM_NAME + "bam.bai"
    threads:
        config['threads']
    shell:
        "samtools index -@ {threads} {input} {output}"
