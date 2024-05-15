import pandas as pd
import numpy as np
import datetime
import re
from snakemake.shell import shell
from snakemake.utils import validate, min_version
from snakemake.io import load_configfile

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Set minimum snakemake version =====                                   │
#└─────────────────────────────────────────────────────────────────────────────┘
min_version("5.7.0")

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Load and process config files =====                                   │
#└─────────────────────────────────────────────────────────────────────────────┘
configfile: "config.yaml"
validate(config, schema="workflow/schema/config.schema.yaml")
config_extra = load_configfile("config_extra.yaml")
# validate(config, schema="workflow/schema/config_extra.schema.yaml")

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Load shared constants and functions =====                             │
#└─────────────────────────────────────────────────────────────────────────────┘
# contains directory structure
include: "workflow/rules/folder_structure.smk"
# contains input functions for rules
# contains other misc functions rules
include: "workflow/rules/input_functions.smk"
include: "workflow/rules/misc_functions.smk"
# contains pipleline rule defitions
# TODO automated parsing for readme?
include: "workflow/rules/pipelines.smk"

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Process metadata =====                                                │
#└─────────────────────────────────────────────────────────────────────────────┘
rnaseq_ext = ["fq", "fastq", "fa", "fasta"]
rnaseq_ext_regex = "|\.".join([e + "$" for e in rnaseq_ext])
paired_read_strings_regex = "|".join(
    [e + "$" for e in config["paired_read_strings"]]
)

Metadata = pd.read_table(config["metadata"], dtype=str)
validate(Metadata, schema="workflow/schema/metadata.schema.yaml")

Metadata["paired"] = np.where(
    Metadata.duplicated(subset=['sample', 'lane'], keep=False), 1, 0
)
Metadata["filename_sans_gzip"] = Metadata["fq"].str.replace(
    ".gz", "", regex=False)
Metadata["filename_sans_ext"] = Metadata["filename_sans_gzip"].str.replace(
    "\." + rnaseq_ext_regex, "", regex=True
)
Metadata["filename_sans_read"] = Metadata["filename_sans_ext"].str.replace(
    paired_read_strings_regex, "", regex=True
)
Metadata = Metadata.sort_values(by=['sample', 'lane', 'read'])

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Global variables =====                                                │
#└─────────────────────────────────────────────────────────────────────────────┘
Samples = list(Metadata["sample"].unique())
Fqs = list(Metadata["fq"].unique())

# string constants based on metadata and config files
DIFFEXP_ANALYSIS = "{}_{}/".format(
    pipelines[config['pipeline']][-1], config["diffexp"]["outdir"])
NOW = str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Load pipleline rules =====                                            │
#└─────────────────────────────────────────────────────────────────────────────┘
# makes sure that certain files that might be needed for analysis are created
# - fasta index file
include: "workflow/rules/preprocess.smk"

# download sra files, skips if sra files are not defined in
include: "workflow/rules/sra.smk"

# run fastqc quality control - load always
include: "workflow/rules/fastqc.smk"

# check whether to load trimming rule and set the fastq input dir for alignment
# accordingly - load conditionally
# TODO maybe simplify the condition for loading the rule
if (
        config['trim']['adapters_single'] != "" or 
        config['trim']['adapters_single'] != "" or 
        config['trim']['quality'] != "" or 
        config['trim']['extra'] != ""
    ):
    include: "workflow/rules/trim.smk"
    FASTQ_INPUT_DIR = TRIMMED_DIR
    workflow_trimming = "yes"
else:
    include: "workflow/rules/skip_trim.smk"
    FASTQ_INPUT_DIR = FASTQ_DIR
    workflow_trimming = "no"

# load pipline rules for selected pipeline from config.yaml 
for rule in pipelines[config['pipeline']]:
    include: RULES_DIR + rule + ".smk"

# TODO fix bam_index rule
# if config['pipeline'] not in ["download_only"]:
#     include: "workflow/rules/bam_index.smk"

# load rules for calculating coverage - load conditionally
if config['coverage']['calculate'] == "yes":
    include: "workflow/rules/coverage.smk"
else:
    include: "workflow/rules/skip_coverage.smk"

# multiqc - load always
include: "workflow/rules/multiqc.smk"

# load rule for creating result archive with processed data - load conditionally
if config['result_archive'] == "yes":
    include: "workflow/rules/result_archive.smk"
else:
    include: "workflow/rules/skip_result_archive.smk"

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Execute before starting the workflow =====                            │
#└─────────────────────────────────────────────────────────────────────────────┘
onstart:
    print()
    print_header("Workflow starting")
    print("Number of fastq files: " + str(len(Fqs)))
    print("Number of samples: " + str(len(Samples)))
    print()
    print_header("Workflow parameters")
    print("Selected pipeline: " + config['pipeline'])
    print("Pipeline steps:")
    for rule in pipelines[config['pipeline']]:
        print(" - " + rule)
    print("Trimming: " + workflow_trimming)
    print("Coverage calculation: " + config['coverage']['calculate'])
    print("Creating result archive: " + config['result_archive'])
    
#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Top level snakemake rule =====                                        │
#└─────────────────────────────────────────────────────────────────────────────┘
rule all:
    input:
        ### check if genomic fasta index file exists
        ### needed for rsem and cufflinks
        config['fasta'] + ".fai",
        ### check if fastq files are present otherwise download sra
        ancient(expand(FASTQ_DIR + "{file}", file=Metadata.fq)),
        ### trim adapters and quality
        get_trim_pe_output_files,
        get_trim_se_output_files,
        get_trim_log_files,
        ### fastqc output files
        get_fastqc_output_files,
        ### align output files
        get_align_output_files,
        get_align_log_files,
        get_bam_index_files,
        ### coverage
        get_coverage_files,
        ### count output files
        get_count_output_files,
        get_count_log_files,
        ### diffexp output files
        get_diffexp_output_files,
        ### multiqc output files
        ### multiqc has problems finding some files
        get_multiqc_output_files,
        ### result archive compressed output
        get_result_archive_output_files


#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Execute after the workflow =====                                      │
#└─────────────────────────────────────────────────────────────────────────────┘
onsuccess:
    from datetime import datetime
    print_header(str(datetime.today().strftime('%Y-%m-%d %H:%M:%S')) + 
    " - RNA-Seq workflow finished successfully!")

onerror:
    from datetime import datetime
    print_header(str(datetime.today().strftime('%Y-%m-%d %H:%M:%S')) + 
    " - There was and ERROR running RNA-Seq workflow!")