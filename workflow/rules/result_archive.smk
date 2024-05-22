def get_result_archive_output_files(wildcards):
    return ARCHIVE_DIR  + config["experiment_name"] + ".tar.gz"

rule result_archive:
    input:
        get_trim_log_files,
        get_fastqc_output_files,
        get_align_output_files,
        get_align_log_files,
        get_coverage_files,
        get_count_output_files,
        get_count_log_files,
        get_diffexp_output_files,
        get_multiqc_output_files
    output:
        ARCHIVE_DIR + config["experiment_name"] + ".tar.gz"
    params:
        RESULT_ARCHIVE_DIRS
    run:
        include_dirs = []
        for item in RESULT_ARCHIVE_DIRS:
            if(os.path.exists(item)):
                include_dirs = include_dirs + [item]
        shell(
            "tar "
            "-czf "    
            "{output} "
            "{include_dirs} "    
        )
        
