"""Wrapper for featurecounts"""

from snakemake.shell import shell

# Metadata = "{snakmake.params.metadata}"

# stranded_str = [
#     featurecounts_stranded(
#         Metadata.query("sample == @x").stranded.dropna().unique()[0]
#     )
#     for x in Samples
# ]

# stranded_str = ",".join(stranded_str)

shell(
    "featureCounts "
    "-a {snakemake.input.gtf} "
    "-o {snakemake.output.counts} "
    "-T {snakemake.threads} "
    #"-s {stranded_str} "
    "{snakemake.params.standard} "
    "{snakemake.params.extra} "
    "{snakemake.input.bam} "
    "2> {snakemake.log}"
)