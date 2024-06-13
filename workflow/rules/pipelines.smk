#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== List of available pipelines =====                                     │
#└─────────────────────────────────────────────────────────────────────────────┘
# requires the components to be basenames of smk rules
pipelines = {
    # only download sra archive or zip files
    "download_only" : {"align": "skip_align", "count": "skip_count", "diffexp": "skip_diffexp"},

    # align only
    "star_only" : {"align": "star", "count": "skip_count", "diffexp": "skip_diffexp"},
    "hisat_only" : {"align": "hisat", "count": "skip_count", "diffexp": "skip_diffexp"},
    "tophat_only": {"align": "tophat", "count": "skip_count", "diffexp": "skip_diffexp"},
    "bowtie2_only": {"align": "bowtie2", "count": "skip_count", "diffexp": "skip_diffexp"},
    "salmon_only": {"align": "salmon", "count": "skip_count", "diffexp": "skip_diffexp"},
    "kallisto_only": {"align": "kallisto", "count": "skip_count", "diffexp": "skip_diffexp"},
    
    # intermediate pipelines
    "featurecounts_only": {"align": "skip_align", "count": "featurecounts", "diffexp": "skip_diffexp"},
    "htseq_only": {"align": "skip_align", "count": "htseq", "diffexp": "skip_diffexp"},
    "stringtie_only": {"align": "skip_align", "count": "stringtie", "diffexp": "skip_diffexp"},
    "stringtie_abundance_only": {"align": "skip_align", "count": "stringtie_abundance", "diffexp": "skip_diffexp"},
    "cufflinks_only": {"align": "skip_align", "count": "cufflinks", "diffexp": "skip_diffexp"},
    "cufflinks_abundance_only": {"align": "skip_align", "count": "cufflinks_abundance", "diffexp": "skip_diffexp"},

    "cufflinks_abundance_cuffdiff": {"align": "skip_align", "count": "cufflinks_abundance", "diffexp": "cuffdiff"},


    "deseq_only": {"align": "skip_align", "count": "featurecounts", "diffexp": "deseq"},
    
    # full diffexp pipelines
    "deseq" : {"align": "star", "count": "featurecounts", "diffexp": "deseq"},
    "cufflinks_cuffdiff": {"align": "skip_align", "count": "cufflinks", "diffexp": "cuffdiff"},
    # "deseq_alt" : {"align": "hisat", "count": "featurecounts", "diffexp": "deseq"},
    # "edger" : {"align": "star", "count": "featurecounts", "diffexp": "edger"},
    # "edger_alt" : {"align": "hisat", "count": "featurecounts", "diffexp": "edger"},
    # "limma" : {"align": "star", "count": "featurecounts", "diffexp": "limma"},
    # "limma_alt" : {"align": "hisat", "count": "featurecounts", "diffexp": "limma"},
    # "stringtie" : {"align": "hisat", "count": "stringtie", "diffexp": "ballgown"},
    # "kallisto" : {"align": "kallisto", "diffexp": "sleuth"},
    # "cuffdiff" : {"align": "tophat", "count": "cufflinks", "diffexp": "cuffdiff"},
    # "cuffdiff_denovo": {"align": "star", "count": "cufflinks_denovo", "diffexp": "cuffdiff"}
}

pipelines_str = ""

for k, v in pipelines.items():
    a =  " - `" + k + "` - "
    b = ", ".join(v.values()) + "\n"
    pipelines_str += a + b

with open("workflow/docs/pipelines.md", "w") as f:
    f.write(pipelines_str)
