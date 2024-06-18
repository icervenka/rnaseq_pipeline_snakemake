"""Wrapper for kallisto alignment."""

from wrapper_functions import arrange_fq_for_align
from snakemake.shell import shell

input_arr = arrange_fq_for_align(
    snakemake.wildcards.sample,
    snakemake.params.metadata,
    snakemake.params.fastq_dir)

if len(input_arr) == 1:
    input_str = "--single " + " ".join(input_arr.iloc[0])
    fragment_info = snakemake.params.fragment_info
elif len(input_arr) == 2:
    input_str =" ".join([ " ".join(x) for x in zip(input_arr.iloc[0], input_arr.iloc[1])])
    fragment_info = ""
else:
    raise ValueError("Too many read types.")

shell(
    "kallisto quant "
    "{snakemake.params.stranded} "
    "{snakemake.params.extra} "
    "{fragment_info} "
    "-b 100 "
    "-t {snakemake.threads} "
    "--genomebam "
    "-g {snakemake.input.gtf} "
    "-i {snakemake.params.index} "
    "-o {snakemake.params.outdir} "
    "{input_str} "
)
