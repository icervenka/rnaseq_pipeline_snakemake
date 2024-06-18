"""Wrapper for fq subsample."""

from snakemake.shell import shell

i1 = snakemake.input[0]
i2 = snakemake.input[1]
o1 = snakemake.output[0]
o2 = snakemake.output[1]

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