def get_cuffnorm_to_matrix(wildcards):
    return rules.counts_to_matrix.output


rule counts_to_matrix:
    input:
        count_table=opj(CUFFNORM_OUTDIR, CUFFNORM_COUNT_TABLE_FILE),
        attr_table=opj(CUFFNORM_OUTDIR, CUFFNORM_ATTR_TABLE_FILE),
        samples_table=opj(CUFFNORM_OUTDIR, CUFFNORM_SAMPLE_TABLE_FILE)
    output:
        opj(COUNT_OUTDIR, COMMON_COUNT_FILE)
    params:
        samples=Samples
    threads:
        1
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/cuffnorm_to_matrix.R"