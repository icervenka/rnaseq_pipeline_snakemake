rule result_archive:
    input:
        get_align_log_files,
        get_count_output_files,
        get_count_log_files,
        get_diffexp_output_files
    output:
        # "archive/" + NOW + "_" + \
        #     config["experiment_name"] + "_result_archive.tar.gz"
        "archive/" + NOW + "_" + config["experiment_name"] + "_result_archive.tar.gz"
    params:
        "diffexp/" + DIFFEXP_OUTDIR
    shell:
        "tar -czf {output} counts {params} logs -C archive"
