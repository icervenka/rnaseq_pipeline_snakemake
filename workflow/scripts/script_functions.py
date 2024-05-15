import os

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    from shutil import which
    return which(name)

##### functions for manipulating/filtering filenames and extensions #####
def parse_filename(string):
    path = os.path.dirname(string)
    filename = os.path.basename(string)
    split1 = os.path.splitext(filename)
    split2 = os.path.splitext(split1[0])
    if(split2[1] == ".fastq"):
        return(split2[0], split1[1][1:])
    else:
        return(split1[0], split1[1][1:])

def extract_name(parsed_filename):
    return(parsed_filename[0])

def extract_extension(parsed_filename):
    return(parsed_filename[1])

##### functions for retrieving fastq file names  #####
def arrange_fq_for_align(sample, metadata, fastq_dir):
    metadata['fq_full'] = fastq_dir + metadata['fq']
    fq_meta = metadata[metadata['sample'] == sample]
    fq_meta = fq_meta.sort_values(by = ["sample", "lane", "read"])
    input_arr = fq_meta.groupby(['read'])['fq_full'].apply(lambda x: x.to_list())
    return input_arr

##### functions for parsing parameters from various tools  #####

def read_command(filename):
    if(filename.endswith(".bz2")):
        return("bunzip2 -c")
    elif(filename.endswith(".gz")):
        return("zcat")
    else:
        return("cat")