# list of available pipelines
# requires the components to be basenames of smk rules
# last item in the list is always script for differential gene expression
pipelines = {
    "cuffdiff" : ['tophat', 'cufflinks', 'cuffdiff'],
    "cuffdiff-denovo: ['star', 'cufflinks-denovo', 'cuffdiff']"
    "stringtie" : ['hisat', 'stringtie', 'ballgown'],
    "deseq" : ['star', 'featurecounts', 'deseq'],
    "deseq-alt" : ['hisat', 'featurecounts', 'deseq'],
    "edger" : ['star', 'featurecounts', 'edger'],
    "edger-alt" : ['hisat', 'featurecounts', 'edger'],
    "limma" : ['star', 'featurecounts', 'limma'],
    "limma-alt" : ['hisat', 'featurecounts', 'limma'],
    "kallisto" : ['kallisto', 'sleuth'],
    "only_download_sra" : ['skip_align', 'skip_count', 'skip_diffexp'],
    "download_align" : ['star', 'skip_count', 'skip_diffexp']
}
