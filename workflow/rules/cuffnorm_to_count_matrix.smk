rule counts_to_matrix:
    input:
        ""
    output:
        ""
    params:
        ""
    script:
        ""

# to make raw read count matrix I need to convert back the fpm size factors
# rule counts_to_matrix:
#     input:
#         expand(rules.count.output.counts, sample=Samples),
#     output:
#         opj(COUNT_OUTDIR, COMMON_COUNT_NAME),
#     conda:
#         CONDA_R_GENERAL_ENV
#     script:
#         "../scripts/featurecounts_gather.R"
