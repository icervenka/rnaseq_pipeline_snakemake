"""Wrapper for MultiQC."""

import os
from snakemake.shell import shell

input_dirs = []

for item in snakemake.input:
    input_dirs = input_dirs + [ os.path.dirname(item) ]

input_dirs = " ".join(list(set(input_dirs)))

shell(
    "multiqc "
    "-f "
    "-d "
    "-dd 1 "
    "-o {snakemake.params.outdir} "
    "-n {snakemake.params.name} "
    "{input_dirs} "
)