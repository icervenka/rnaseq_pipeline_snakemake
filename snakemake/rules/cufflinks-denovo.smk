rule count:
    input:
        bam = rules.rename_bam.output,
        gtf = config["gtf"],
        fasta = config["fasta"]
    output:
        counts = COUNT_OUTDIR + "{sample}/transcripts.gtf"
    params:
        extra = cufflinks_params
    threads:
        config["threads"]
    run:
        shell(
            "cufflinks "
            "-g {input.gtf} "
            "--multi-read-correct "
            "--upper-quartile-norm "
            "--library-type {params.extra} "
            "-o {COUNT_OUTDIR}{sample} "
            "{input.bam} "
        )
