"""Wrapper for star alignment."""

import os
from script_functions import read_command, arrange_fq_for_align
from snakemake.shell import shell

input_arr = arrange_fq_for_align(
    snakemake.wildcards.sample,
    snakemake.params.metadata,
    snakemake.params.fastq_dir)

input_arr = [','.join(x) for x in input_arr.to_list()]
input_str = " ".join(input_arr)

outprefix = os.path.dirname(snakemake.output.bam) + "/"
read_cmd = read_command(input_arr[0])

shell(
    "STAR "
    "{snakemake.params.extra} "
    "--runThreadN {snakemake.threads} "
    "--genomeDir {snakemake.params.index} "
    "--readFilesIn {input_str} "
    "--outFileNamePrefix {outprefix} "
    "--readFilesCommand {read_cmd} "
    "--outSAMstrandField intronMotif "
    "--outSAMtype BAM SortedByCoordinate "
    "--outFilterType BySJout "
    "--limitBAMsortRAM 20000000000 "
)
