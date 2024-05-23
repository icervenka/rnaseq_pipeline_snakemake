# TODO verify the ballgown filename

def get_count_output_files(wildcards):
    # return expand(COUNT_OUTDIR + "{sample}/transcripts.gtf", sample=Samples) +
    #     expand(COUNT_OUTDIR + "{sample}/gene_counts.tsv", sample=Samples) +
    #     expand(COUNT_OUTDIR + "{sample}/{ballgown}.ctab", sample=Samples,
    #         ballgown=BALLGOWN_INPUT_FILES) +
    #     [COUNT_OUTDIR + "counts.tsv",
    #         COUNT_OUTDIR + "fpkm.tsv",
    #         COUNT_OUTDIR + "samples_combined.tsv"]
    # return expand(rules.count.output.counts, sample=Samples)
    return (
        [rules.counts_to_matrix.output.counts_gene,
        rules.counts_to_matrix.output.counts_transcript, 
        rules.gather_data.output.tpm,
        rules.gather_data.output.fpkm,
        rules.gather_data.output.samples_combined] +
        expand(rules.count.output.gtf, sample=Samples) + 
        expand(rules.count.output.gtf, sample=Samples)
    )

def get_count_log_files(wildcards):
    return []

rule count:
    input:
        bam=rules.align_out.output,
        gtf=config["gtf"]
    output:
        counts=COUNT_OUTDIR + "{sample}/gene_counts.txt",
        gtf=COUNT_OUTDIR + "{sample}/{sample}.gtf",
        ballgown=expand(COUNT_OUTDIR + "{{sample}}/{file}", file=BALLGOWN_INPUT_FILES)
    params:
        #stranded=stringtie_stranded,
        #standard=stringtie_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    threads:
        config["threads"]
    conda:
        CONDA_COUNT_GENERAL_ENV
    shell:
        """
        stringtie \
        {params.extra} \
        -p {threads} \
        -G {input.gtf} \
        -B \
        -e \
        -A {output.counts}\
        -o {output.gtf} \
        {input.bam} 
        """
# -e \
# "{params.stranded} "

rule counts_to_matrix:
    input:
        expand(rules.count.output.counts, sample=Samples)
    output:
        counts_gene=COUNT_OUTDIR + "counts.csv",
        counts_transcript=COUNT_OUTDIR + "transcript_counts.csv"
    params:
        input_dir=COUNT_OUTDIR,
        length=75
    threads:
        1
    conda:
        CONDA_COUNT_GENERAL_ENV
    shell:
        """
        python3 workflow/scripts/prepDE.py \
        -l {params.length} \
        -i {params.input_dir} \
        -g {output.counts_gene} \
        -t {output.counts_transcript} \
        """

rule gather_data:
    input:
        expand(rules.count.output.counts, sample=Samples)
    output:
        tpm=COUNT_OUTDIR + "tpm.txt",
        fpkm=COUNT_OUTDIR + "fpkm.txt",
        samples_combined=COUNT_OUTDIR + "samples_combined.txt"
    conda:
        CONDA_R_GENERAL_ENV
    script:
        "../scripts/stringtie_count_gather.R"
