from snakemake.shell import shell

inp = snakemake.params[0][0]

shell(
    "while true; "
    "do "
        "if ! vdb-validate {snakemake.output}; "
        "then "
            " prefetch -o {snakemake.output} {inp} || true ; "
        "else "
            "break; "
        "fi "
    "done "
)