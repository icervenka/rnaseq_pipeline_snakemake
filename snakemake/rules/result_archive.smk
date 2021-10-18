
# TODO change input to something more dynamic and "complete"
rule result_archive:
    input:
        LOG_DIR + DIFFEXP_ANALYSIS + "diffexp.completed"
    output:
        # "archive/" + NOW + "_" + \
        #     config["experiment_name"] + "_result_archive.tar.gz"
        "archive/" + config["experiment_name"] + "_result_archive.tar.gz"
    params:
        "diffexp/" + DIFFEXP_OUTDIR
    shell:
        "tar -czf {output} counts {params} logs -C archive"
