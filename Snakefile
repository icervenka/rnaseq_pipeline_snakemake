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
include: "snakemake/rules/functions.smk"
include: "snakemake/rules/pipelines.smk"

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

##### string constants based on metadata and config files #####
#SAMPLES_BASENAME = [extract_name(parse_filename(x)) for x in Samples]
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

include: "snakemake/rules/bam_index.smk"

if config['coverage']['calculate'] == "yes":
    include: "snakemake/rules/coverage.smk"
else:
    include: "snakemake/rules/skip_coverage.smk"

#include: "snakemake/rules/multiqc.smk"
include: "snakemake/rules/result_archive.smk"

##### top level snakemake rule #####
rule all:
    input:
        # check if genomic fasta index file exists
        # needed for rsem and cufflinks
        config['fasta'] + ".fai",
        # check if fastq files are present otherwise download sra
        expand(FASTQ_DIR + "{file}", file=Metadata.fq),
        # get_sra_download_files,
        # fastqc output files
        #expand(LOG_DIR + "qc/{file}.html", file=Metadata.fq),
        #expand(LOG_DIR + "qc/{file}_fastqc.zip", file=Metadata.fq),
        # align output files
        get_align_output_files,
        get_align_log_files,
        # get_bam_index_files,
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
        "archive/" + NOW + "_" + config["experiment_name"] + "_result_archive.tar.gz"
