import pandas as pd
import datetime
import re
from snakemake.shell import shell
from snakemake.utils import validate, min_version
from snakemake.io import load_configfile

##### set minimum snakemake version #####
min_version("5.7.0")

##### load config and sample sheets #####
# pipelines, string constants and functions
include: "snakemake/rules/common.smk"
include: "snakemake/rules/input_functions.smk"
include: "snakemake/rules/pipelines.smk"

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

# config
configfile: "config.yaml"
validate(config, schema="snakemake/schema/config.schema.yaml")
config_extra = load_configfile("config_extra.yaml")

# metadata
Metadata = pd.read_table(config["metadata"], dtype=str)
validate(Metadata, schema="snakemake/schema/metadata.schema.yaml")
Metadata = Metadata.sort_values(by=['sample', 'lane', 'read'])

Samples = list(Metadata["sample"].unique())
Fqs = list(Metadata["fq"].unique())
Basenames = [ extract_name(parse_filename(x)) for x in Fqs ]
print(Basenames)

##### string constants based on metadata and config files #####
DIFFEXP_ANALYSIS = "{}_{}/".format(
    pipelines[config['pipeline']][-1], config["diffexp"]["outdir"])
NOW = str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))

##### load pipleline rules #####
include: "snakemake/rules/preprocess.smk"
include: "snakemake/rules/sra.smk"
#include: "snakemake/rules/fastqc.smk"
include: "snakemake/rules/trim.smk"

for rule in pipelines[config['pipeline']]:
    include: RULES_DIR + rule + ".smk"

# if config['pipeline'] not in ["download_only"]:
#     include: "snakemake/rules/bam_index.smk"

if config['coverage']['calculate'] == "yes":
    include: "snakemake/rules/coverage.smk"
else:
    include: "snakemake/rules/skip_coverage.smk"

#include: "snakemake/rules/multiqc.smk"

# TODO archive requires all directories to be present
# include: "snakemake/rules/result_archive.smk"

##### top level snakemake rule #####
rule all:
    input:
        # check if genomic fasta index file exists
        # needed for rsem and cufflinks
        config['fasta'] + ".fai",
        # check if fastq files are present otherwise download sra
        expand(FASTQ_DIR + "{file}", file=Metadata.fq),
        get_trim_output_files,
        # fastqc output files
        #expand(LOG_DIR + "qc/{file}.html", file=Metadata.fq),
        #expand(LOG_DIR + "qc/{file}_fastqc.zip", file=Metadata.fq),
        # align output files
        get_align_output_files,
        get_align_log_files,
        #get_bam_index_files,
        # coverage get_align_output_files
        get_coverage_files,
        # count output files
        get_count_output_files,
        get_count_log_files,
        # diffexp output files
        get_diffexp_output_files,
        # multiqc output files
        # multiqc has problems finding some files
        #LOG_DIR + "qc/" + re.sub("\s+", "", config["experiment_name"]) + ".html",
        # result archive compressed output
        #"archive/" + NOW + "_" + config["experiment_name"] + "_result_archive.tar.gz"
