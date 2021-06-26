import os

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

def read_command(filename):
    if(filename.endswith(".bz2")):
        return("bunzip2 -c")
    elif(filename.endswith(".gz")):
        return("zcat")
    else:
        return("cat")

def get_fq(wildcards):
    fq = metadata.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    return fq
