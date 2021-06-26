from snakemake.shell import shell

fq1 = snakemake.input.sample[0]
assert fq1 is not None, "input-> fq1 is a required input parameter"
fq1 = (
    [snakemake.input.sample[0]]
    if isinstance(snakemake.input.sample[0], str)
    else snakemake.input.sample[0]
)

fq2 = snakemake.input.sample[1] if len(snakemake.input.sample) > 1 else None
if fq2:
    fq2 = (
        [snakemake.input.sample[1]]
        if isinstance(snakemake.input.sample[1], str)
        else snakemake.input.sample[1]
    )
    assert len(fq1) == len(
        fq2
    ), "input-> equal number of files required for fq1 and fq2"

input_str_fq1 = ",".join(fq1)
if fq is None:
    input_str = "-U " + input_str_fq1
else:
    input_str = "-1 " + input_str_fq1 + " -2 " + ",".join(fq2)

# TODO check if it works with gz and bz2 files
shell(
    "hisat2 "
    "{snakemake.params.extra} "
    "-p {snakemake.threads} "
    "-q "
    "-x {snakemake.input.index} "
    "{input_str} "
    "-S {snakemake.output.bam} "
    ">> {snakemake.output.log} 2>&1 "
)
