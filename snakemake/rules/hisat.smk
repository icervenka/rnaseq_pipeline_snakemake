rule align:
    input:
        sample = get_fq,
        index = config["index"]
    output:
        sam = temp(ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".sam"),
        log = ALIGN_LOG_OUTDIR + "{sample}/hisat.log"
    params:
        extra = config['align']['hisat_extra']
    threads:
        config["threads"]
    script:
        "../scripts/hisat_wrapper.py"

rule sort_sam:
    input:
        rules.align.sam
    output:
        ALIGN_OUTDIR + "{sample}/" + COMMON_BAM_NAME + ".bam",
    params:
        compression = 9
        flag = "0x2"
    threads:
        config['threads']
    run:
        shell(
            "samtools view "
            "-uSh "
            #    "-f {params.flag} "
            "{input} "
            "| "
            "samtools sort "
            "-l {params.compression} "
            "-@ {threads} "
            "--output-fmt BAM "
            "-o {output} "
        )

rule move_log:
    input:
        rules.align.output.log
    output:
        ALIGN_LOG_OUTDIR + "{sample}/hisat.log"
    shell:
        "mv {input} {output}"

rule all_align:
    input:
        expand(ALIGN_OUTDIR + "{sample}/" +
               COMMON_BAM_NAME + ".bam", sample=Samples),
        expand(ALIGN_LOG_OUTDIR + "{sample}/hisat.log", sample=Samples)
    output:
        touch(ALIGN_LOG_OUTDIR + "align.completed")
