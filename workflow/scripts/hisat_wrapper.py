"""Wrapper for hisat2 alignment."""

from script_functions import arrange_fq_for_align
from snakemake.shell import shell

input_arr = arrange_fq_for_align(
    snakemake.wildcards.sample,
    snakemake.params.metadata,
    snakemake.params.fastq_dir)

input_arr = [','.join(x) for x in input_arr.to_list()]

if len(input_arr) == 1:
    input_str = "-U " + input_arr[0]
elif len(input_arr) == 2:
    input_str = "-1 " + input_arr[0] + " -2 " + input_arr[1]
else:
    raise ValueError("Too many read types.")

# TODO check if it works with gz and bz2 files
shell(
    "hisat2 "
    "{snakemake.params.extra} "
    "-p {snakemake.threads} "
    "-q "
    "-x {snakemake.params.index} "
    "{input_str} "
    "-S {snakemake.output.sam} "
    "--novel-splicesite-outfile {snakemake.output.splicesite} "
    "--met-file {snakemake.output.met} "
    ">> {snakemake.output.log} 2>&1 "
)
