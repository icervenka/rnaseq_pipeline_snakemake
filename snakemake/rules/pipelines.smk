# list of available pipelines
# requires the components to be basenames of smk rules
# last item in the list is always script for differential gene expression
pipelines = {
    "cuffdiff" : ['tophat', 'cufflinks', 'cuffdiff'],
    "cuffdiff_denovo: ['star', 'cufflinks-denovo', 'cuffdiff']"
    "stringtie" : ['hisat', 'stringtie', 'ballgown'],
    "deseq" : ['star', 'featurecounts', 'deseq'],
    "deseq_alt" : ['hisat', 'featurecounts', 'deseq'],
    "edger" : ['star', 'featurecounts', 'edger'],
    "edger_alt" : ['hisat', 'featurecounts', 'edger'],
    "limma" : ['star', 'featurecounts', 'limma'],
    "limma_alt" : ['hisat', 'featurecounts', 'limma'],
    "kallisto" : ['kallisto', 'sleuth'],
    "download_only" : ['skip_align', 'skip_count', 'skip_diffexp'],
    "star_only" : ['star', 'skip_count', 'skip_diffexp'],
    "featurecounts_only": ['skip_align', 'featurecounts', 'skip_diffexp'],
    "deseq_only": ['skip_align', 'featurecounts', 'deseq'],
}
