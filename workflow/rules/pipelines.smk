#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== List of available pipelines =====                                     │
#└─────────────────────────────────────────────────────────────────────────────┘
# requires the components to be basenames of smk rules
pipelines = {
    # only download sra archive or zip files
    "download_only" : {"align": "skip_align", "count": "skip_count", "diffexp": "skip_diffexp"},

    # align only pipelines
    "star_only" : {"align": "star", "count": "skip_count", "diffexp": "skip_diffexp"},
    "hisat_only" : {"align": "hisat", "count": "skip_count", "diffexp": "skip_diffexp"},
    "tophat_only": {"align": "tophat", "count": "skip_count", "diffexp": "skip_diffexp"},
    # "bowtie2_only": {"align": "bowtie2", "count": "skip_count", "diffexp": "skip_diffexp"},
    "salmon_only": {"align": "salmon", "count": "skip_count", "diffexp": "skip_diffexp"},
    "kallisto_only": {"align": "kallisto", "count": "skip_count", "diffexp": "skip_diffexp"},
    
    # count only pipelines
    "featurecounts_only": {"align": "skip_align", "count": "featurecounts", "diffexp": "skip_diffexp"},
    "htseq_only": {"align": "skip_align", "count": "htseq", "diffexp": "skip_diffexp"},
    "stringtie_only": {"align": "skip_align", "count": "stringtie", "diffexp": "skip_diffexp"},
    "cufflinks_only": {"align": "skip_align", "count": "cufflinks", "diffexp": "skip_diffexp"},

    # diffexp only pipelines
    "deseq_only": {"align": "skip_align", "count": "skip_count", "diffexp": "deseq"},
    "cuffdiff_only": {"align": "skip_align", "count": "skip_count", "diffexp": "deseq"},

    # intermediate pipelines
    "featurecounts_deseq": {"align": "skip_align", "count": "featurecounts", "diffexp": "deseq"},
    "htseq_deseq": {"align": "skip_align", "count": "htseq", "diffexp": "deseq"},
    "stringtie_deseq": {"align": "skip_align", "count": "stringtie", "diffexp": "deseq"},
    
    "cufflinks_cuffdiff": {"align": "skip_align", "count": "cufflinks", "diffexp": "cuffdiff"},
    "stringtie_cuffdiff": {"align": "skip_align", "count": "stringtie", "diffexp": "cuffdiff"},
    
    # full pipelines
    "deseq" : {"align": "star", "count": "featurecounts", "diffexp": "deseq"},
    "deseq" : {"align": "hisat", "count": "htseq", "diffexp": "deseq"},
    "edger" : {"align": "star", "count": "featurecounts", "diffexp": "edger"},
    "edger_alt" : {"align": "hisat", "count": "htseq", "diffexp": "edger"},
    "limma" : {"align": "star", "count": "featurecounts", "diffexp": "limma"},
    "limma_alt" : {"align": "hisat", "count": "htseq", "diffexp": "limma"},
    "stringtie" : {"align": "hisat", "count": "stringtie", "diffexp": "ballgown"},
    "kallisto" : {"align": "kallisto", "diffexp": "sleuth"},
    "cufflinks": {"align": "tophat", "count": "cufflinks", "diffexp": "cuffdiff"},
}

pipelines_str = ""

for k, v in pipelines.items():
    a =  " - `" + k + "` - "
    b = ", ".join(v.values()) + "\n"
    pipelines_str += a + b

with open("workflow/docs/pipelines.md", "w") as f:
    f.write(pipelines_str)
