def get_count_output_files(wildcards):
    return expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples) +
        [ COUNT_OUTDIR + "merged.gtf" ]

def get_count_log_files(wildcards):
    return COUNT_OUTDIR + "gtf_assembly.txt"

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
            "-G {input.gtf} "
            "--frag-bias-correct {input.fasta} "
            "{params.extra} "
            "-o {COUNT_OUTDIR}{sample} "
            "{input.bam} "
        )

rule assemble:
    input:
        expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples)
    output:
        COUNT_OUTDIR + "gtf_assembly.txt"
    shell:
        "printf '%s\n' {input} >> {output}"

rule merge:
    input:
        merged=rules.assemble.output,
        gtf=config['gtf']
    output:
        COUNT_OUTDIR + "merged.gtf"
    params:
        extra=cuffmerge_params
    threads:
        config['threads']
    run:
        shell(
            "export PS1=; "
            "source /usr/local/bin/miniconda3/etc/profile.d/conda.sh; "
            "conda activate tophat; "
            "cuffmerge "
            "-p {threads} "
            "-g {input.gtf} "
            "{params.extra} "
            "-o {COUNT_OUTDIR} "
            "{input.merged} "
        )
