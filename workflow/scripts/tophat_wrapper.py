from snakemake.shell import shell

input_arr = arrange_fq_for_align(
    snakemake.input.sample,
    snakemake.params.metadata,
    snakemake.params.fastq_dir)

input_arr = [','.join(x) for x in input_arr.to_list()]
input_str = " ".join(input_arr)

shell(
    "export PS1=; "
    "source /usr/local/bin/miniconda3/etc/profile.d/conda.sh; "
    "conda activate tophat; "
    "tophat2 "
    "{snakemake.params.extra} "
    "-p {snakemake.threads} "
    "-G {snakemake.input.gtf} "
    "-o {snakemake.params.align_outdir} "
    "{snakemake.input.index} "
    "{input_str} "
)
