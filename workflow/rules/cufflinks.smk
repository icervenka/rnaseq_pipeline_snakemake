# TODO add cuffnorm

def get_count_output_files(wildcards):
    return expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples) +
        [ COUNT_OUTDIR + "merged.gtf" ]

def get_count_log_files(wildcards):
    return [ COUNT_LOG_OUTDIR + "gtf_assembly.txt" ]

rule assemble:
    input:
        bam=rules.align_out.output,
        gtf=config["gtf"],
        fasta=config["fasta"]
    output:
        COUNT_OUTDIR + "{sample}/transcripts.gtf"
    params:
        extra=cufflinks_params
    threads:
        config["threads"]
    run:
        if config_extra['count']['cufflinks_mode'] == "classic":
            mode="-G {input.gtf} "
        elif config_extra['count']['cufflinks_mode'] == "guided":
            mode="-g {input.gtf} "
        elif config_extra['count']['cufflinks_mode'] == "denovo":
            mode=""
        else:
            mode="-G {input.gtf} "

        shell(
            "cufflinks "
            "{mode} "
            "--frag-bias-correct {input.fasta} "
            "{params.extra} "
            "-o {COUNT_OUTDIR}{sample} "
            "{input.bam} "
        )

rule gather_gtf:
    input:
        expand(rules.assemble.output, sample=Samples)
    output:
        COUNT_LOG_OUTDIR + "gtf_assembly.txt"
    shell:
        "printf '%s\n' {input} >> {output}"

# TODO check if gtf should be specified also for denovo assembly
rule merge:
    input:
        gathered=rules.gather_gtf.output,
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
            "{input.gathered} "
        )
