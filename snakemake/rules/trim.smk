ruleorder: cutadapt_pe > cutadapt_se

rule cutadapt_se:
    input:
        get_single_end_cutadapt
    output:
        TRIMMED_DIR + "{fq}.fastq"
    params:
        adapters=config['trim']['adapters_single'],
        quality=config['trim']['quality'],
        extra=config['trim']['extra']
    threads: 4
    log:
        LOG_DIR + "trim/{fq}.txt"
    run:
        if params.adapters == "" and params.quality == "" and params.extra == "":
            shell(
                "ln -s {input} {output.fastq}; "
                "printf 'No trimming performed\n' > {output.qc}"
            )
        else:
            shell(
                "cutadapt "
                "{params.adapters} "
                "{params.quality} "
                "{params.extra} "
                "-j {threads} "
                "-o {output.fastq} "
                "{input} "
                "> {log} "
            )

rule cutadapt_pe:
    input:
        get_paired_end_cutadapt
    output:
        expand(TRIMMED_DIR + "{{fq}}{read}.fastq", read=['_1', "_2"])
    params:
        adapters=config['trim']['adapters_paired'],
        quality=config['trim']['quality'],
        extra=config['trim']['extra']
    threads: 4
    log:
        LOG_DIR + "trim/{fq}.txt"
    run:
        if params.adapters == "" and params.quality == "" and params.extra == "":
            shell(
                "ln -s {input[0]} {output[0]}; "
                "ln -s {input[1]} {output[1]}; "
                "printf 'No trimming performed\n' > {output.qc}"
            )
        else:
            shell(
                "cutadapt "
                "{params.adapters} "
                "{params.quality} "
                "{params.extra} "
                "-j {threads} "
                "-o {output[0]} "
                "-p {output[1]} "
                "{input[0]} "
                "{input[1]} "
                "> {output.qc} "
            )
