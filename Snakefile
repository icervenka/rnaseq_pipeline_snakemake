import pandas as pd
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

##### top level snakemake rule #####
rule all:
    input:
        #LOG_DIR + "sra.completed",
        expand(FASTQ_DIR + "{file}", file=Metadata.fq),
        LOG_DIR + "qc.completed",
        LOG_DIR + "align.completed",
        LOG_DIR + "count.completed",
        LOG_DIR + DIFFEXP_ANALYSIS + "diffexp.completed"

##### load remaining pipleline rules #####
include: "snakemake/rules/sra.smk"

# TODO
if config['trim'] == "yes":
    include: "snakemake/rules/trim.smk"

# TODO
if config['coverage'] == "yes":
    include: "snakemake/rules/coverage.smk"

for rule in pipelines[config['pipeline']]:
    include: RULES_DIR + rule + ".smk"

if config['result_archive'] == "yes":
    include: "snakemake/rules/result_archive.smk"

include: "snakemake/rules/qc.smk"
