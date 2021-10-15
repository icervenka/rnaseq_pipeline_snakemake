# TODO reflect new folder structure
rule result_archive:
    input:
        DESEQ2_OUTDIR + "reports/report.html",
        "logs/" + config["experiment_name"] + ".html"
    output:
        "archive/" + NOW + "_" + \
            config["experiment_name"] + "_result_archive.tar.gz"
    shell:
        "tar -czf {output} counts diffexp logs metadata.tsv config.yaml Snakefile -C archive"
