ruleorder: cutadapt_pe > cutadapt_se

# TODO fix extensions in trimming

def get_trim_pe_output_files(wildcards):
    return expand(TRIMMED_DIR + "{filename}{read}" + ".fastq.gz", 
        filename=list(
            Metadata[Metadata["paired"] == 1]["filename_sans_read"].unique()
        ), 
        read=config["paired_read_strings"]
    )

def get_trim_se_output_files(wildcards):
    return expand(TRIMMED_DIR + "{filename}.fastq.gz", 
        filename=list(
            Metadata[Metadata["paired"] == 0]["filename_sans_read"].unique()
        )
    )

def get_trim_log_files(wildcards):
    return expand(LOG_DIR + "trim/{filename}.txt", 
        filename=list(
            Metadata["filename_sans_read"].unique()
        )
    )


rule cutadapt_se:
    input:
        FASTQ_DIR + "{filename}.fastq.gz"
    output:
        TRIMMED_DIR + "{filename}.fastq.gz"
    params:
        adapters=config['trim']['adapters_single'],
        quality=config['trim']['quality'],
        extra=config['trim']['extra']
    threads: 4
    log:
        LOG_DIR + "trim/{filename}.txt"
    run:
        if params.adapters == "" and params.quality == "" and params.extra == "":
            shell(
                "ln {input} {output}; "
                "printf 'No trimming performed\n' > {log}"
            )
        else:
            shell(
                "cutadapt "
                "{params.adapters} "
                "{params.quality} "
                "{params.extra} "
                "-j {threads} "
                "-o {output} "
                "{input} "
                "> {log} "
            )

rule cutadapt_pe:
    input:
        expand(FASTQ_DIR + "{{filename}}{read}" + ".fastq.gz", 
            read=config["paired_read_strings"]
        )
    output:
        expand(TRIMMED_DIR + "{{filename}}{read}" + ".fastq.gz", 
            read=config["paired_read_strings"]
        )
    params:
        adapters=config['trim']['adapters_paired'],
        quality=config['trim']['quality'],
        extra=config['trim']['extra']
    threads: 4
    log:
        LOG_DIR + "trim/{filename}.txt"
    run:
        i1 = input[0]
        i2 = input[1]
        o1 = output[0]
        o2 = output[1]
        print(wildcards)
        if params.adapters == "" and params.quality == "" and params.extra == "":
            shell(
                "ln {i1} {o1}; "
                "ln {i2} {o2}; "
                "printf 'No trimming performed\n' > {log}"
            )
        else:
            shell(
                "cutadapt "
                "{params.adapters} "
                "{params.quality} "
                "{params.extra} "
                "-j {threads} "
                "-o {o1} "
                "-p {o2} "
                "{i1} "
                "{i2} "
                "> {log} "
            )
