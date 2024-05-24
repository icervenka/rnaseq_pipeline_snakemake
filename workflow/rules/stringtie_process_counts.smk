 def get_stringtie_processed_output_files(wildcards):
    return (
        [rules.counts_to_matrix.output.counts_gene,
        rules.counts_to_matrix.output.counts_transcript, 
        rules.gather_data.output.tpm,
        rules.gather_data.output.fpkm,
        rules.gather_data.output.samples_combined] 
    )

rule counts_to_matrix:
    input:
        expand(rules.count.output.counts, sample=Samples)
    output:
        counts_gene=COUNT_OUTDIR + COMMON_COUNT_NAME,
        counts_transcript=COUNT_OUTDIR + COMMON_TRANSCRIPT_COUNT_NAME
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
        tpm=COUNT_OUTDIR + STRINGTIE_TPM_FILE,
        fpkm=COUNT_OUTDIR + STRINGTIE_FPKM_FILE,
        samples_combined=COUNT_OUTDIR + STRINGTIE_COMBINED_FILE
    params:
        samples=expand("{sample}", sample=Samples)
    conda:
        CONDA_R_GENERAL_ENV
    script:
        "../scripts/stringtie_count_gather.R"