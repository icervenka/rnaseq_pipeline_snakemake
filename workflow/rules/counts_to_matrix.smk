rule counts_to_matrix:
    input:
        files=get_tximport_files,
        gtf=config["gtf"]
    output:
        #opj(COUNT_OUTDIR, COMMON_COUNT_NAME)
        "counts/count_matrix.txt"
    params:
        metadata=config["metadata"],
        samples=Samples,
        species=config["species"],
        pipeline=pipeline.values(),
        tool=config_extra["tximport_count_matrix"]["tool"],
        ids_in=config["diffexp"]["input_gene_ids"],
        out_type="counts",
        extra=config_extra["tximport_count_matrix"]["extra"]
    threads:
        1
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/tximport.R"