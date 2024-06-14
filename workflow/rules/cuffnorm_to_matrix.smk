# TODO Do the files need to be encoded in _tools_filenames?
rule counts_to_matrix:
    input:
        count_table=opj(CUFFNORM_OUTDIR, "genes.count_table"),
        attr_table=opj(CUFFNORM_OUTDIR, "genes.attr_table"),
        samples_table=opj(CUFFNORM_OUTDIR, "samples.table.")
    output:
        opj(CUFFNORM_OUTDIR, COMMON_COUNT_NAME)
    params:
        samples=Samples
    threads:
        1
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/cuffnorm_to_matrix.R"