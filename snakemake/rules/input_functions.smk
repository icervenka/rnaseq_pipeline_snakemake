def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    from shutil import which
    return which(name)

##### functions for retrieving and dumping sra files #####
def get_sra(wildcards):
    sra = Metadata.query('sra == @wildcards.sra').sra.dropna().unique()
    return sra

def get_single_end_sra(wildcards):
    m = Metadata[~Metadata.duplicated(subset=['sra'], keep=False).astype(bool)]
    m = m.query('sra == @wildcards.sra').sra.unique()
    return SRA_DIR + m

def get_paired_end_sra(wildcards):
    m = Metadata[Metadata.duplicated(subset=['sra']).astype(bool)]
    m = m.query('sra == @wildcards.sra').sra.unique()
    return SRA_DIR + m

##### functions for getting names of fastq files for trimming #####
# def get_single_end_cutadapt(wildcards):
#     m = Metadata[Metadata.paired == 0]
#     m = m.query('filename_sans_read == @wildcards.filename')
#     return list(FASTQ_INPUT_DIR + m.fq)

# def get_paired_end_cutadapt(wildcards):
#     m = Metadata[Metadata.paired == 1]
#     m = m.query('filename_sans_read == @wildcards.filename')
#     print(m)
#     return list(FASTQ_INPUT_DIR + m.fq)

##### functions for retrieving fastq file names  #####
def get_fq(wildcards):
    m = Metadata.query('sample == @wildcards.sample').dropna()
    return FASTQ_INPUT_DIR + m.fq

def arrange_fq_for_align(samples, metadata, fastq_dir):
    metadata['fq_full'] = fastq_dir + metadata['fq']
    fq_meta = metadata[metadata['fq_full'].isin(samples)]
    input_arr = fq_meta.groupby(['read'])['fq_full'].apply(lambda x: x.to_list())
    return input_arr

##### functions for parsing parameters from various tools  #####
def featurecounts_stranded(stranded):
    def stranded_switch(x):
        select = {
            "no": "0",
            "yes": "1",
            "reverse": "2"
        }
        return(select.get(x, "0"))
    return stranded_switch(stranded)

def featurecounts_params(wildcards):
    param_string = ""
    param_string += "-M " if config["count"]["multimap"] == "yes" else ""
    param_string += "-O " if config["count"]["overlap"] == "yes" else ""
    param_string += config_extra["count"]["featurecounts_extra"] + " "
    return(param_string)

def htseq_params(wildcards):
    s = Metadata.query('sample == @wildcards.sample').stranded.dropna().unique()[0]
    param_string = ""
    param_string += s + " "
    if config["count"]["overlap"] == "yes":
        param_string += "--nonunique=all "
        if config["count"]["multimap"] == "no":
            param_string += "--secondary-alignments=ignore "
    param_string += config_extra["count"]["htseq_extra"] + " "
    return(param_string)

def kallisto_params(wildcards):
    def stranded_switch(x):
        select = {
            "no": " ",
            "yes": "--fr-stranded ",
            "reverse": "--rf-stranded "
        }
        return(select.get(x, " "))

    s = Metadata.query('sample == @wildcards.sample').stranded.dropna().unique()
    param_string = ""
    param_string += stranded_switch(s)
    param_string += config_extra['align']['kallisto_extra']
    return(param_string)

def stringtie_stranded(wildcards):
    select = {
        "no": "",
        "yes": "--rf ",
        "reverse": "--fr "
    }
    s = Metadata.query('sample == @wildcards.sample').stranded.dropna().unique()[0]
    return(select.get(s, ""))

# TODO change numbers for yes/no/reverse
def cufflinks_params(wildcards):
    def stranded_switch(x):
        select = {
            "1": "ff-firststrand ",
            "2": "ff-secondstrand ",
            "3": "ff-unstranded ",
            "4": "fr-firststrand ",
            "5": "fr-secondstrand ",
        }
        return(select.get(x, "ff-firststrand"))

    s = Metadata.query('sample == @wildcards.sample').stranded.dropna().unique()[0]

    param_string = ""
    param_string += "--library-type" + stranded_switch(s)
    param_string += config_extra['count']['cufflinks_extra']
    return(param_string)

def cuffmerge_params(wildcards):
    param_string = ""
    param_string += config_extra['count']['cuffmerge_extra']
    return(param_string)

##### input function for creating cuffdiff data payloads  #####

def get_cuffdiff_data():
    cond = config['diffexp']['design']
    cond = cond.replace('~', '')
    cond = cond.strip()

    ref_levels = config['diffexp']['ref_levels']

    if len(ref_levels) > 0:
        ref_levels = "".join(ref_levels.split())
        groups = ref_levels.split(',')
        labels_str = ref_levels
    else:
        groups = Metadata[cond].unique()
        labels_str = ','.join(groups)

    files = []
    for group in groups:
        sample_files = Metadata[Metadata[cond].eq(group)]['sample'].to_list()
        files.append(','.join([ALIGN_OUTDIR + x + "/" +
                               COMMON_BAM_NAME + ".bam" for x in sample_files]))

    files_str = " ".join(files)

    return (labels_str, files_str)
