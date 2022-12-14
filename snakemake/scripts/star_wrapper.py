from functions import read_command, arrange_fq_for_align
import os
from snakemake.shell import shell

input_arr = arrange_fq_for_align(
    snakemake.input.sample,
    snakemake.params.metadata,
    snakemake.params.fastq_dir)

# meta = snakemake.params.metadata
# meta['fq_full'] = snakemake.params.fastq_dir + meta['fq']
# fq_meta = snakemake.params.metadata[snakemake.params.metadata['fq_full'].isin(snakemake.input.sample)]
# input_arr = fq_meta.groupby(['read'])['fq_full'].apply(lambda x: x.to_list())
input_arr = [','.join(x) for x in input_arr.to_list()]
input_str = " ".join(input_arr)

outprefix = os.path.dirname(snakemake.output.bam) + "/"
print(input_arr)
read_cmd = read_command(input_arr[0])

shell(
    "STAR "
    "{snakemake.params.extra} "
    "--runThreadN {snakemake.threads} "
    "--genomeDir {snakemake.input.index} "
    "--readFilesIn {input_str} "
    "--outFileNamePrefix {outprefix} "
    "--readFilesCommand {read_cmd} "
    "--outSAMtype BAM SortedByCoordinate "
    "--outFilterType BySJout "
    "--limitBAMsortRAM 20000000000 "
)
