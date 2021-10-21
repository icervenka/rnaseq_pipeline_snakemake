from snakemake.shell import shell

# fq1 = snakemake.input.sample[0]
# assert fq1 is not None, "input-> fq1 is a required input parameter"
# fq1 = (
#     [snakemake.input.sample[0]]
#     if isinstance(snakemake.input.sample[0], str)
#     else snakemake.input.sample[0]
# )
#
# fq2 = snakemake.input.sample[1] if len(snakemake.input.sample) > 1 else None
# if fq2:
#     fq2 = (
#         [snakemake.input.sample[1]]
#         if isinstance(snakemake.input.sample[1], str)
#         else snakemake.input.sample[1]
#     )
#     assert len(fq1) == len(
#         fq2
#     ), "input-> equal number of files required for fq1 and fq2"
#
# input_str_fq1 = ",".join(fq1)
# if fq is None:
#     input_str = "-U " + input_str_fq1
# else:
#     input_str = "-1 " + input_str_fq1 + " -2 " + ",".join(fq2)

meta = snakemake.params.metadata
meta['fq_full'] = snakemake.params.fastq_dir + meta['fq']
fq_meta = snakemake.params.metadata[snakemake.params.metadata['fq_full'].isin(snakemake.input.sample)]
input_arr = fq_meta.groupby(['read'])['fq_full'].apply(lambda x: x.to_list())
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
    "-S {snakemake.output.bam} "
    ">> {snakemake.output.log} 2>&1 "
)
