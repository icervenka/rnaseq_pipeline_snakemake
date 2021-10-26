from snakemake.shell import shell

# # TODO: redirect to log
# meta = Metadata
# meta['fq_full'] = FASTQ_DIR + meta['fq']
# fq_meta = meta[meta['fq_full'].isin(input.sample)]
# print(fq_meta)
# # TODO: test if ordering is ok on multilane fastqfiles
# input_arr = fq_meta.sort_values(by=['sample']).groupby(['read'])['fq_full'].apply(lambda x: x.to_list())
# input_arr = [' '.join(x) for x in input_arr.to_list()]
# input_str = " ".join(input_arr)
# print(input_str)

input_arr = arrange_fq_for_align(
    snakemake.input.sample,
    snakemake.params.metadata,
    snakemake.params.fastq_dir)

if len(input_arr) == 1:
    input_str = "--single " + " ".join(input_arr[0])
    fragment_info = snakemake.params.single_extra
elif len(input_arr) == 2:
    input_str = " ".join(zip(input_arr[0], input_arr[1]))
    fragment_info = ""
else:
    raise ValueError("Too many read types.")


# if(fq_meta['read'].unique().size == 1):
#     mode = "--single "
#     fragment_info = "-l " + str(params.length_mean) + " -s " + str(params.length_variance)
# else:
#     mode = " "
#     fragment_info = " "

# stranded = kallisto_stranded(fq_meta['stranded'].iloc[0])
# print(stranded)
shell(
    "kallisto quant "
    "-b 100 "
    "-t {snakemake.threads} "
    "--genomebam "
    "-g {snakemake.input.gtf} "
    "-i {snakemake.input.index} "
    "-o {ALIGN_OUTDIR}{snakemake.wildcards.sample} "
    "{snakemake.params.extra} "
    "{input_str} "
)
