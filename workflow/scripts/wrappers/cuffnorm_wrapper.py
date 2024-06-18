"""Wrapper for cuffnorm."""

from snakemake.shell import shell
from os import path
import warnings

Metadata = snakemake.params.metadata
samples_stranded = Metadata["stranded"].unique().tolist()

if len(samples_stranded) > 1:
    warnings.warn("Samples have differing strandedness specified, " + 
                  "defaulting to unstranded library type.")
    stranded = "fr-unstranded"
else:
    stranded = samples_stranded[0]

shell(
    """
    cuffnorm \
    -q \
    --no-update-check \
    --output-format simple-table \
    --library-norm-method geometric \
    -p {snakemake.threads} \
    -o {snakemake.params.outdir} \
    -L {snakemake.params.labels} \
    --library-type {stranded} \
    {snakemake.params.extra} \
    {snakemake.input.gtf} \
    {snakemake.input.bam} \
    > {snakemake.log} 2>&1; 
    rm {snakemake.params.outdir}/run.info
    """
)
