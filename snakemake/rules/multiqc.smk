import re

def get_multiqc_output_files(wildcards):
    return(LOG_DIR +
           "multiqc/" +
           re.sub("\s+", "", config["experiment_name"]) +
           ".html")

rule multiqc:
    input:
        get_trim_log_files,
        get_fastqc_output_files,
        get_align_log_files,
        get_count_log_files
    output:
        html = LOG_DIR + "multiqc/" + \
            re.sub("\s+", "", config["experiment_name"]) + ".html",
    params:
        name = re.sub("\s+", "", config["experiment_name"]),
        outdir = LOG_DIR + "multiqc/"
    run:
        input_dirs = []
        for item in input:
            input_dirs = input_dirs + [ os.path.dirname(item) ]
        input_dirs = " ".join(list(set(input_dirs)))
        print(input_dirs)
        shell(
            "multiqc "
            "-f "
            "-d "
            "-dd 1 "
            "-o {params.outdir} "
            "-n {params.name} "
            "{input_dirs} "
        )
        
