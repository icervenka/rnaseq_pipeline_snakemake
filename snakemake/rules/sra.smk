ruleorder: fastq_dump_paired > fastq_dump_single
# uses prefetch-validate combo because fasterq-dump has no way of
# checking if the file was downloaded correctly due to network errors

rule sra_download:
    input:
        config['metadata']
    output:
        temp(SRA_DIR + "{sra}")
    params:
        get_sra
    threads:
        1
    run:
        inp = params[0][0]
        shell(
            "while true; "
            "do "
                "if ! vdb-validate {output}; "
                "then "
                    " prefetch -o {output} {inp} || true ; "
                "else "
                    "break; "
                "fi "
            "done "
        )

rule fastq_dump_single:
    input:
        get_single_end_sra
    output:
        FASTQ_DIR + "{sra}.fastq"
    threads:
        config['threads']
    run:
        shell(
            "fasterq-dump "
            "-e {threads} "
            "-s "
            "{input} "
            "-o {output}"
        )

rule fastq_dump_paired:
    input:
        get_paired_end_sra
    output:
        expand(FASTQ_DIR + "{{sra}}{read}.fastq", read=['_1', "_2"])
    threads:
        config['threads']
    run:
        shell(
            "fasterq-dump "
            "-e {threads} "
            "-3 "
            "{input} "
            "-o " + FASTQ_DIR + "{wildcards.sra}.fastq"
        )

rule compress_fastq:
    input:
        FASTQ_DIR + "{file}.fastq"
    output:
        FASTQ_DIR + "{file}.fastq.gz"
    threads:
        config['threads']
    run:
        has_pigz = is_tool("pigz")
        program = "pigz" if has_pigz is not None else "gzip"
        run_threads = ("-p " + str(threads)) if has_pigz else ""
        shell(
            "{program} "
            "--best "
            "{run_threads} "
            "{input} "
        )
