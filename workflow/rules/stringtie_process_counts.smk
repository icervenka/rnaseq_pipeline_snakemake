def get_stringtie_processed_output_files(wildcards):
    return [rules.stringtie_gather.output.tpm,
        rules.stringtie_gather.output.fpkm,
        rules.stringtie_gather.output.samples_combined]


rule stringtie_gather:
    input:
        expand(rules.count.output.counts, sample=Samples)
    output:
        tpm=opj(STRINGTIE_OUTDIR, STRINGTIE_TPM_FILE),
        fpkm=opj(STRINGTIE_OUTDIR, STRINGTIE_FPKM_FILE),
        samples_combined=opj(STRINGTIE_OUTDIRR, STRINGTIE_COMBINED_FILE)
    params:
        samples=expand("{sample}", sample=Samples)
    conda:
        CONDA_R_GENERAL_ENV
    script:
        "../scripts/stringtie_gather.R"


# TODO compare with tximport
##  remove last row with colsums
# rule counts_to_matrix:
#     input:
#         expand(rules.count.output.counts, sample=Samples)
#     output:
#         counts_gene=opj(COUNT_OUTDIR, COMMON_COUNT_FILE,)
#         counts_transcript=opj(COUNT_OUTDIR, COMMON_TRANSCRIPT_COUNT_FILE)
#     params:
#         input_dir=COUNT_OUTDIR,
#         length=75
#     threads:
#         1
#     conda:
#         CONDA_COUNT_GENERAL_ENV
#     shell:
#         """
#         python3 workflow/scripts/prepDE.py \
#         -l {params.length} \
#         -i {params.input_dir} \
#         -g {output.counts_gene} \
#         -t {output.counts_transcript} \
#         """