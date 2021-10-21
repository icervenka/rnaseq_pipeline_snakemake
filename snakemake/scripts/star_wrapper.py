import os
from snakemake.shell import shell

# def is_tool(name):
#     """Check whether `name` is on PATH and marked as executable."""
#
#     from shutil import which
#     return which(name)

# def read_command(filename):
#     if(filename.endswith(".bz2")):
#         return("bunzip2 -c")
#     elif(filename.endswith(".gz")):
#         if(is_tool("pigz") is not None):
#             return(f"pigz -d -p {snakemake.threads} -c")
#         else:
#             return("zcat")
#     else:
#         return("cat")

meta = snakemake.params.metadata
meta['fq_full'] = snakemake.params.fastq_dir + meta['fq']
fq_meta = snakemake.params.metadata[snakemake.params.metadata['fq_full'].isin(snakemake.input.sample)]
input_arr = fq_meta.groupby(['read'])['fq_full'].apply(lambda x: x.to_list())
input_arr = [','.join(x) for x in input_arr.to_list()]
input_str = " ".join(input_arr)

outprefix = os.path.dirname(snakemake.output.bam) + "/"
read_cmd = read_command(fq_meta["fq"].iloc[0])

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
