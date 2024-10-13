def get_coverage_files(wildcards):
    config['coverage']['split_strands'] == "no":
        return expand(COVERAGE_OUTDIR + "{sample}.bedgraph", sample=Samples)
    elif config['coverage']['split_strands'] == "yes":
        return expand(COVERAGE_OUTDIR + "{sample}_{strand}.bedgraph",
            sample=Samples, strand=['forward', 'reverse'])

if config['coverage']['split_strands'] == "no"
    rule bam_coverage:
        input:
            rules.align_out.output[0]
        output:
            COVERAGE_OUTDIR + "{sample}.bedgraph"
        params:
            normalize_using=config['coverage']['normalize_using'],
            bin_size=config['coverage']['bin_size'],
        threads:
            config['threads']
        run:
            shell(
                "bamCoverage "
                "-p {threads} "
                "-bs {params.bin_size} "
                "-of bedgraph "
                "--normalizeUsing {params.normalize_using} "
                "-b {input} "
                "-o {output} "
            )
elif config['coverage']['split_strands'] == "yes":
    rule bam_coverage:
        input:
            rules.align_out.output[0]
        output:
            forward=COVERAGE_OUTDIR + "{sample}_forward.bedgraph",
            reverse=COVERAGE_OUTDIR + "{sample}_reverse.bedgraph"
        params:
            normalize_using=config['coverage']['normalize_using'],
            bin_size=config['coverage']['bin_size'],
        threads:
            config['threads']
        run:
            shell(
                "bamCoverage "
                "-p {threads} "
                "-bs {params.bin_size} "
                "-of bedgraph "
                "--normalizeUsing {params.normalize_using} "
                "--filterRNAstrand forward "
                "-b {input} "
                "-o {output.forward} "
            )
            shell(
                "bamCoverage "
                "-p {threads} "
                "-bs {params.bin_size} "
                "-of bedgraph "
                "--normalizeUsing {params.normalize_using} "
                "--filterRNAstrand reverse "
                "-b {input} "
                "-o {output.reverse} "
            )
else:
    raise ValueError("Only yes or no allowed as value for splitting strands.")

# rule invert_minus_bedgraph:
#     input:
#         rules.bam_coverage.output.minus
#     output:
#         COVERAGE_OUTDIR + "{sample}_minus_inverted.bedgraph
#     shell:
#         "awk ‘$4 *= -1’ {input} > {output}"

rule bedgraph_to_tdf:
    input:
        COVERAGE_OUTDIR + "{bedgraph}.bedgraph"
    output:
        COVERAGE_OUTDIR + "{bedgraph}.tdf"
    params:
        config['fasta']
    run:
        shell(
            "igvtools toTDF "
            "-z 7 "
            "{input} "
            "{output} "
            "{params.fasta} "
        )
