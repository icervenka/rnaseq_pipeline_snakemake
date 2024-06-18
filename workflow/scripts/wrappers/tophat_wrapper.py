"""Wrapper for tophat alignment."""

from wrapper_functions import arrange_fq_for_align
from snakemake.shell import shell

input_arr = arrange_fq_for_align(
    snakemake.wildcards.sample,
    snakemake.params.metadata,
    snakemake.params.fastq_dir)

input_arr = [','.join(x) for x in input_arr.to_list()]
input_str = " ".join(input_arr)

shell(
    "tophat2 "
    "--library-type {snakemake.params.stranded} "
    "{snakemake.params.extra} "
    "-p {snakemake.threads} "
    "-G {snakemake.input.gtf} "
    "-o {snakemake.params.align_outdir} "
    "{snakemake.params.index} "
    "{input_str} "
)
