"""Wrapper for cuffdiff."""

from snakemake.shell import shell
from os import path
import warnings

# create a modify tilda function
Metadata = snakemake.params.metadata
Samples = snakemake.params.samples
samples_stranded = Metadata["stranded"].unique().tolist()

group = snakemake.params.diffexp["design"].replace('~', '').strip()
ref_levels = snakemake.params.diffexp["reference_levels"][group]
fdr = snakemake.params.diffexp["fdr"]
_ds = snakemake.params.ds
_tf = snakemake.params.tf
extra = snakemake.params.diffexp["extra"]


if len(samples_stranded) > 1:
    warnings.warn("Samples have differing strandedness specified, " + 
                  "defaulting to unstranded library type.")
    stranded = "fr-unstranded"
else:
    stranded = samples_stranded[0]
    
if len(ref_levels) == 0:
    labels = Metadata[group].unique()
elif len(ref_levels) == 1:
    labels = Metadata[group].unique()
    labels = [ref_levels] + [ x for x in labels if ref_levels != x ]
else:
    labels = ref_levels

files = []
for item in labels:
    sample_files = Metadata[Metadata[group].eq(item)]['sample'].to_list()
    files.append(
        ','.join(
            [ path.join(_ds["ALIGN_OUTDIR"], x, _tf["COMMON_BAM_NAME"] + ".bam") for x in sample_files ]
        )
    )

files_str = " ".join(files)
labels_str = ",".join(labels)

shell(
    """
    cuffdiff \
    --no-update-check \
    -q \
    {extra} \
    -p {snakemake.threads} \
    -o {snakemake.params.outdir} \
    --frag-bias-correct {snakemake.params.fasta} \
    --library-type {stranded} \
    --FDR {fdr} \
    -u {snakemake.input.gtf} \
    --labels {labels_str} \
    {files_str} \
    > {snakemake.log} 2>&1
    """
)
