"""Wrapper for fq subsample."""

from snakemake.shell import shell

i1 = snakemake.input[0]
i2 = snakemake.input[1]
o1 = snakemake.output[0]
o2 = snakemake.output[1]

shell(
    """
    "fq subsample \
    {params.proportion} \
    --r1-dst {o1} \
    --r2-dst {o2} \
    {i1} \
    {i2} \
    "> {log} 2>&1 "
    """
)