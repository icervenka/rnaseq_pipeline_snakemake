#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== List of available pipelines =====                                     │
#└─────────────────────────────────────────────────────────────────────────────┘
# requires the components to be basenames of smk rules
# last item in the list is always script for differential gene expression
pipelines = {
    # only download sra archive or zip files
    "download_only" : ['skip_align', 'skip_count', 'skip_diffexp'],

    # align only
    "star_only" : ['star', 'skip_count', 'skip_diffexp'],
    "hisat_only" : ['hisat', 'skip_count', 'skip_diffexp'],
    "tophat_only": ['tophat', 'skip_count', 'skip_diffexp'],
    "bowtie2_only": ['bowtie2', 'skip_count', 'skip_diffexp'],
    "salmon_only": ['salmon', 'skip_count', 'skip_diffexp'],
    "kallisto_only": ['kallisto', 'skip_count', 'skip_diffexp'],
    
    # intermediate pipelines
    "featurecounts_only": ['skip_align', 'featurecounts', 'skip_diffexp'],
    "htseq_only": ['skip_align', 'htseq', 'skip_diffexp'],
    "stringtie_only": ['skip_align', 'stringtie', 'skip_diffexp'],
    "stringtie_expression_only": ['skip_align', 'stringtie_expression', 'skip_diffexp'],
    "cufflinks_only": ['skip_align', 'cufflinks', 'skip_diffexp'],

    "deseq_only": ['skip_align', 'featurecounts', 'deseq'],
    
    # full diffexp pipelines
    "deseq" : ['star', 'featurecounts', 'deseq']
    # "deseq_alt" : ['hisat', 'featurecounts', 'deseq'],
    # "edger" : ['star', 'featurecounts', 'edger'],
    # "edger_alt" : ['hisat', 'featurecounts', 'edger'],
    # "limma" : ['star', 'featurecounts', 'limma'],
    # "limma_alt" : ['hisat', 'featurecounts', 'limma'],
    # "stringtie" : ['hisat', 'stringtie', 'ballgown'],
    # "kallisto" : ['kallisto', 'sleuth'],
    # "cuffdiff" : ['tophat', 'cufflinks', 'cuffdiff'],
    # "cuffdiff_denovo": ['star', 'cufflinks_denovo', 'cuffdiff']
}

pipelines_str = ""

for k, v in pipelines.items():
    # print(" - `" + k + "` - ")
    # print(", ".join(v) + "\n")
    a =  " - `" + k + "` - "
    b = ", ".join(v) + "\n"
    pipelines_str += a + b

with open("workflow/docs/pipelines.md", "w") as f:
    f.write(pipelines_str)
