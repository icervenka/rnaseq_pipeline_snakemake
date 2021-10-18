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

Samples = list(Metadata["sample"].unique())
Fqs = list(Metadata["fq"].unique())

##### string constants based on metadata and config files #####
#SAMPLES_BASENAME = [extract_name(parse_filename(x)) for x in Samples]
DIFFEXP_ANALYSIS = "{}_{}/".format(
    pipelines[config['pipeline']][-1], config["diffexp"]["outdir"])
NOW = str(datetime.datetime.now())

print(expand(FASTQ_DIR + "{file}", file=Metadata.fq))
##### top level snakemake rule #####
rule all:
    input:
        #LOG_DIR + "sra.completed",
        expand(FASTQ_DIR + "{file}", file=Metadata.fq),
        # for fastqc rule
        expand(LOG_DIR + "qc/{file}.html", file=Metadata.fq),
        expand(LOG_DIR + "qc/{file}_fastqc.zip", file=Metadata.fq),
        # # for multiqc rule
        # # multiqc has problems finding some files
        # LOG_DIR + "qc/" + re.sub("\s+", "", config["experiment_name"]) + ".html",
        # #LOG_DIR + "align.completed",
        expand(ALIGN_OUTDIR+"{sample}/"+COMMON_BAM_NAME+".bam", sample=Samples),
        # # TODO this will only work for star logs, package other logs into function probably
        expand(ALIGN_LOG_OUTDIR+"{sample}/{log}", sample=Samples, log=STAR_LOGFILES),
        LOG_DIR + "count.completed",
        LOG_DIR + DIFFEXP_ANALYSIS + "diffexp.completed",
        # # for result archive
        "archive/" + config["experiment_name"] + "_result_archive.tar.gz"

##### load remaining pipleline rules #####
include: "snakemake/rules/sra.smk"
include: "snakemake/rules/fastqc.smk"

# TODO
if config['trim'] == "yes":
    include: "snakemake/rules/trim.smk"

# TODO
if config['coverage'] == "yes":
    include: "snakemake/rules/coverage.smk"

for rule in pipelines[config['pipeline']]:
    include: RULES_DIR + rule + ".smk"

#include: "snakemake/rules/multiqc.smk"
include: "snakemake/rules/result_archive.smk"

# if config['result_archive'] == "yes":
#     include: "snakemake/rules/result_archive.smk"
