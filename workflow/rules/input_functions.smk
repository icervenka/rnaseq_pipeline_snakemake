def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    from shutil import which
    return which(name)

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Functions for retrieving and dumping sra files =====                  │
#└─────────────────────────────────────────────────────────────────────────────┘
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

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Functions for retrieving fastq file names =====                       │
#└─────────────────────────────────────────────────────────────────────────────┘
def get_fq(wildcards):
    m = Metadata.query('sample == @wildcards.sample').dropna()
    return FASTQ_CURRENT_DIR + m.fq


def get_single_fq(wildcards, fastq_dir):
    return f"{fastq_dir}{wildcards.filename}{wildcards.ext}"


def get_paired_fq(wildcards, fastq_dir):
    return expand(f"{{fastq_dir}}{{wildcards.filename}}{read}{{ext}}", 
            read=config["paired_read_strings"])


#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Functions for acessing metadata columns =====                         │
#└─────────────────────────────────────────────────────────────────────────────┘
def get_metadata(what, paired = 1):
    if paired != 1 and paired != 0:
        raise ValueError("Value for 'paired' argument can only be 0 or 1.")
    return Metadata[Metadata["paired"] == paired][what].unique()


#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Functions for parsing parameters from various tools =====             │
#└─────────────────────────────────────────────────────────────────────────────┘
# TODO rsem has overlapping values with featurecounts, will lead to errors
# tools with empty string in unstranded protocol that also need an argument name
# need to have the passing of the argument name handled by respecive stranded 
# function instead of the  wrapper. Currently: hisat
def get_strandedness(stranded, outtool):
    d = {
        "featurecounts": ["0", "1", "2"],
        "htseq": ["yes", "no", "reverse"],
        "cufflinks": ["fr-unstranded", "fr-secondstrand", "fr-firststrand"],
        "stringtie": ["", "--fr", "--rf"],
        # "rsem": ["0.5", "1", "0"],
        "hisat_se": ["", "F", "R"], 
        "hisat_pe": ["", "FR", "RF"],
        "tophat": ["fr-unstranded", "fr-secondstrand", "fr-firststrand"],
        "kallisto": ["", "--fr-stranded", "--rf-stranded"],
        "salmon_se": ["U", "SF", "SR"],
        "salmon_pe": ["IU", "ISF", "ISR"]
    }
    df = pd.DataFrame(d)
    mask = df.map(lambda x: stranded in str(x)).any(axis=1)
    res = df[mask]

    try:
        res = res[outtool].tolist()[0]
    except IndexError:
        print("No value for strandedness returned, is it " + 
            "specified correctly in the metadata file")
        
    return res

def get_sample_strandedness(wildcards):
    s = Metadata.query('sample == @wildcards.sample').stranded.dropna().unique()[0]
    return s

def stranded_param(wildcards, tool):
    return get_strandedness(get_sample_strandedness(wildcards), tool) + " "

# featurecounts ----------------------------------------------------------------
def featurecounts_multi_stranded(wildcards):
    s = Metadata.query('sample == @wildcards.sample').stranded.dropna()
    stranded_arr = [ get_strandedness(x, "featurecounts") for x in s ]
    stranded_str = ",".join(stranded_arr)
    return "-s " + stranded_str + " "


def featurecounts_params(wildcards):
    param_string = ""
    param_string += "-M " if config["count"]["multimap"] else ""
    param_string += "-O " if config["count"]["overlap"] else ""
    param_string += config_extra["count"]["featurecounts_extra"] + " "
    return(param_string)

# htseq ------------------------------------------------------------------------
# def htseq_stranded(wildcards):
#     s = Metadata.query('sample == @wildcards.sample').stranded.dropna().unique()[0]
#     return "-s " + get_strandedness(s, "htseq") + " "


def htseq_params(wildcards):
    param_string = ""
    if config["count"]["overlap"] == "yes":
        param_string += "--nonunique=all "
    if not config["count"]["multimap"]:
        param_string += "--secondary-alignments=ignore "
    return(param_string)

# stringtie --------------------------------------------------------------------
# def stringtie_stranded(wildcards):
#     s = Metadata.query('sample == @wildcards.sample').stranded.dropna().unique()[0]
#     return get_strandedness(m, "stringtie") + " "

def stringtie_params(wildcards):
    return ""

# cufflinks --------------------------------------------------------------------
# def cufflinks_stranded(wildcards):
#     s = Metadata.query('sample == @wildcards.sample').stranded.dropna().unique()[0]
#     return "--library-type" + get_strandedness(s, "cufflinks") + " "


def cufflinks_params(wildcards):
    # -u/–multi-read-correct
    # –max-multiread-fraction <0.0-1.0>
    # –no-effective-length-correction
    if config_extra['count']['cufflinks_mode'] == "classic":
        mode="-G {input.gtf} "
    elif config_extra['count']['cufflinks_mode'] == "guided":
        mode="-g {input.gtf} "
    elif config_extra['count']['cufflinks_mode'] == "denovo":
        mode=""
    else:
        mode="-G {input.gtf} "

    return(param_string)


def cuffmerge_params(wildcards):
    param_string = ""
    param_string += config_extra['count']['cuffmerge_extra']
    return(param_string)

# hisat ------------------------------------------------------------------------
def hisat_stranded(wildcards):
    r = Metadata.query('sample == @wildcards.sample').read.dropna().unique()
    outtool = "hisat_pe" if len(r.tolist()) > 1 else "hisat_se"
    return get_strandedness(get_sample_strandedness(wildcards), outtool)

# tophat -----------------------------------------------------------------------
# def tophat_stranded(wildcards):
#     s = Metadata.query('sample == @wildcards.sample').stranded.dropna().unique()[0]
#     return "--library-type" + get_strandedness(s, "tophat") + " "

# kallisto ---------------------------------------------------------------------
# def kallisto_stranded(wildcards):
#     s = Metadata.query('sample == @wildcards.sample').stranded.dropna().unique()[0]
#     return get_strandedness(s, "kallisto") + " "

# salmon -----------------------------------------------------------------------
def salmon_stranded(wildcards):
    r = Metadata.query('sample == @wildcards.sample').read.dropna().unique()
    outtool = "salmon_pe" if len(r.tolist()) > 1 else "salmon_se"
    return get_strandedness(get_sample_strandedness(wildcards), outtool)

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Input function for creating cuffdiff data payloads =====              │
#└─────────────────────────────────────────────────────────────────────────────┘
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


#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Other input/param functions =====                                     │
#└─────────────────────────────────────────────────────────────────────────────┘
def get_subsample_proportion(n):
    if type(n) != int and type(n) != float:
        raise ValueError("Incorrect value for propotion of sampled reads." + 
            "Must be either int or float.")
    else:
        if n <= 1:
            return "-p " + str(n)
        else:
            return "-n " + str(n)

def has_extra_config(conf, conf_extra):
    conf = conf.strip()
    if len(conf) > 0:
        if conf_extra[conf] == None or len(conf_extra[conf]) == 0:
            return conf + " "
            
        else:
            return conf_extra[conf]
    else:
        return " "
