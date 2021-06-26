# TODO needs fasta index to be present: samtools faidx

rule count:
    input:
        bam = rules.rename_bam.output,
        gtf = config["gtf"],
        fasta = config["fasta"]
    output:
        counts = COUNT_OUTDIR + "{sample}/transcripts.gtf"
    threads:
        config["threads"]
    run:
        shell(
            "cufflinks "
            "-G {input.gtf} "
            "--multi-read-correct "
            "--upper-quartile-norm "
            "--frag-bias-correct {input.fasta} "
            "-o {COUNT_OUTDIR}{sample} "
            "{input.bam} "
        )

rule assemble:
    input:
        expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples)
    output:
        COUNT_OUTDIR + MERGE_OUTDIR + "gtf_assembly.txt"
    shell:
        "printf '%s\n' {input} >> {output}"

rule merge:
    input:
        merged = rules.assemble.output,
        gtf = config['gtf']
    output:
        COUNT_OUTDIR + MERGE_OUTDIR + "merged.gtf"
    threads:
        config['threads']
    run:
        shell(
            "export PS1=; "
            "source /usr/local/bin/miniconda3/etc/profile.d/conda.sh; "
            "conda activate tophat; "
            "cuffmerge "
            "-o {COUNT_OUTDIR}{MERGE_OUTDIR} "
            "-p {threads} "
            "-g {input.gtf} "
            "{input.merged} "
        )

rule all_count:
    input:
        expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples),
        COUNT_OUTDIR + MERGE_OUTDIR + "gtf_assembly.txt",
        COUNT_OUTDIR + MERGE_OUTDIR + "merged.gtf"
    output:
        touch(LOG_DIR + "count.completed")
