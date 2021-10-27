def get_count_output_files(wildcards):
    return expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples)

rule count:
    input:
        bam=rules.rename_bam.output,
        gtf=config["gtf"],
        fasta=config["fasta"]
    output:
        counts=COUNT_OUTDIR + "{sample}/transcripts.gtf"
    params:
        extra=cufflinks_params
    threads:
        config["threads"]
    run:
        shell(
            "cufflinks "
            "-g {input.gtf} "
            "{params.extra} "
            "-o {COUNT_OUTDIR}{sample} "
            "{input.bam} "
        )
