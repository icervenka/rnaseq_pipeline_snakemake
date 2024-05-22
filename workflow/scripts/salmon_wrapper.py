from script_functions import arrange_fq_for_align
from snakemake.shell import shell

input_arr = arrange_fq_for_align(
    snakemake.wildcards.sample,
    snakemake.params.metadata,
    snakemake.params.fastq_dir)


# input_arr = [','.join(x) for x in input_arr.to_list()]

if len(input_arr) == 1:
    input_str = "-r " + " ".join(input_arr.iloc[0])
    fragment_info = snakemake.params.fragment_info
elif len(input_arr) == 2:
    input_str ="-1 " + " ".join(input_arr.iloc[0]) + " -2 " + " ".join(input_arr.iloc[1])
    # input_str = " ".join(zip(input_arr.iloc[0], input_arr.iloc[1]))
    fragment_info = ""
else:
    raise ValueError("Too many read types.")

shell(
        "salmon quant "
        "--softclipOverhangs "
        "--validateMappings "
        "{snakemake.params.extra} "
        "{fragment_info} "
        "-i {snakemake.params.index} "
        "-l {snakemake.params.stranded} "
        "{input_str} "
        "-o {snakemake.params.outdir} "
)
