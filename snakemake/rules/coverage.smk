# TODO infer from infer experiment

rule sam_index:
    input:
        rules.rename_bam.output[0]
    output:
        ALIGN_OUTDIR + "{sample}.bai"
    shell:
        "samtools index {input} {output}"

rule bam_coverage:
    input:
        sample = rules.rename_bam.output[0]
        fasta = config['fasta']
    output:
        plus = COVERAGE_OUTDIR + "{sample}_plus.bedgraph
        minus = temp(COVERAGE_OUTDIR + "{sample}_minus.bedgraph)
    threads:
        config['threads']
    run:
        shell(
            "bamCoverage "
            "-p {threads} "
            "-bs 1 "
            "-of bedgraph "
            "--normalizeUsing BPM "
            "--filterRNAstrand reverse "
            "-b {input.sample} "
            "-o {output.plus} "
        )

        shell(
            "bamCoverage "
            "-p {threads} "
            "-bs 1 "
            "-of bedgraph "
            "--normalizeUsing BPM "
            "--filterRNAstrand forward "
            "-b {input.sample} "
            "-o {output.minus} "
        )

rule invert_minus_bedgraph:
    input:
        rules.bam_coverage.output.minus
    output:
        COVERAGE_OUTDIR + "{sample}_minus_inverted.bedgraph
    shell:
        "awk ‘$4 *= -1’ {input} > {output}"

rule bedgraph_to_tdf:
    input:

        rules.bam_coverage.output.plus,
        rules.invert_minus_bedgraph.output[0]
    output:
        COVERAGE_OUTDIR + {}
    run:
        shell(
            "igvtools toTDF "
            "-z 7 "
            ""
            ""
            "{input.fasta} "
        )
