import pandas as pd
import numpy as np
import datetime
import re
from os import path
from snakemake.shell import shell
from snakemake.utils import validate, min_version
from snakemake.io import load_configfile

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Misc setup =====                                                      │
#└─────────────────────────────────────────────────────────────────────────────┘
# Is used a lot in all the rules, this might save some typing
opj = path.join

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
# load global string constants:
# - directory structure for input output and workflow
# - filenames created by different tools used in workflow
# - conda environment names
RULES_DIR = opj("workflow", "rules")
include: opj(RULES_DIR, "_globals.smk")
include: opj(RULES_DIR, "_tool_filenames.smk")
# contains misc functions and input functions for rules
include: opj(RULES_DIR, "functions.smk")
# contains pipleline rule defitions
include: opj(RULES_DIR, "pipelines.smk")

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Process metadata =====                                                │
#└─────────────────────────────────────────────────────────────────────────────┘
# load and validate metadata ---------------------------------------------------
Metadata = pd.read_table(config["metadata"], dtype=str)
validate(Metadata, schema="workflow/schema/metadata.schema.yaml")

# add other useful file columns ------------------------------------------------
# regex to determine compression of files
compression_ext = ["gz", "bz2"]
compression_ext_regex = "(\\." + "|\.".join([e + "$" for e in compression_ext]) + ")"

# regex to determine file extensions of files
rnaseq_ext = ["fq", "fastq", "fa", "fasta"]
rnaseq_ext_regex = "(\\." + "|\.".join([e + "$" for e in rnaseq_ext]) + ")"

# regex to deteremine which reads are paired
paired_read_strings_regex = "|".join(
    [e + "$" for e in config["paired_read_strings"]]
)

Metadata["paired"] = np.where(
    Metadata.duplicated(subset=["sample", "lane"], keep=False), 1, 0
)

# Add additional columns to metadata 
Metadata["filename_sans_compression"] = Metadata["fq"].str.replace(
    compression_ext_regex, "", regex=True)

Metadata["filename_sans_ext"] = Metadata["filename_sans_compression"].str.replace(
    rnaseq_ext_regex, "", regex=True
)

Metadata["filename_sans_read"] = Metadata["filename_sans_ext"].str.replace(
    paired_read_strings_regex, "", regex=True
)

Metadata["filename_compression"] = Metadata["fq"].str.extract(compression_ext_regex)
Metadata["filename_ext"] = Metadata["filename_sans_compression"].str.extract(rnaseq_ext_regex)
Metadata["filename_full_ext"] = Metadata["filename_ext"] + Metadata["filename_compression"]

Metadata = Metadata.sort_values(by=["sample", "lane", "read"])

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Global sample and fastq file variables =====                          │
#└─────────────────────────────────────────────────────────────────────────────┘
Samples = list(Metadata["sample"].unique())
Fqs = list(Metadata["fq"].unique())

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Load pipleline rules =====                                            │
#└─────────────────────────────────────────────────────────────────────────────┘
# download sra files, skips if sra files are not defined
include: opj(RULES_DIR, "sra.smk")

# makes sure that certain files that might be needed for analysis are created
# - fasta index file
include: opj(RULES_DIR, "preprocess.smk")

# run fastqc quality control - load always
include: opj(RULES_DIR, "fastqc.smk")

# load pipline rules for selected pipeline from config.yaml 
pipeline = pipelines[config["pipeline"]]
for key, rule in pipeline.items():
    include: opj(RULES_DIR, rule + ".smk")

# load rules for calculating coverage - load conditionally
if config["coverage"]["calculate"]:
    include: opj(RULES_DIR, "coverage.smk")
else:
    include: opj(RULES_DIR, "skip_coverage.smk")

# multiqc - load always
include: opj(RULES_DIR, "multiqc.smk")

# load rule for creating result archive with processed data - load conditionally
if config["result_archive"]:
    include: opj(RULES_DIR, "result_archive.smk")
else:
    include: opj(RULES_DIR, "skip_result_archive.smk")

#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Execute before starting the workflow =====                            │
#└─────────────────────────────────────────────────────────────────────────────┘
onstart:
    print()
    print_header("Workflow starting")
    print("Number of fastq files:", len(Fqs))
    print("Number of samples:", len(Samples))
    print()
    print_header("Workflow parameters")
    print("Selected pipeline:", config["pipeline"])
    print("Pipeline steps:")
    for key, rule in pipeline.items():
        print(" -",  key , ":", rule)
    
    if is_set_subsample(config["preprocess"]["subsample"]):
        print("Subsampling:", config["preprocess"]["subsample"])
    else:
        print("Subsampling:  False")
    
    if is_set_trimmer(config["trim"]["trimmer"]):
        print("Trimming: Using", config["trim"]["trimmer"])
    else:
        print("Trimming: False")
    
    print("Coverage calculation: " + str(config["coverage"]["calculate"]))
    print("Creating result archive: " + str(config["result_archive"]))
    
#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Top level snakemake rule =====                                        │
#└─────────────────────────────────────────────────────────────────────────────┘
rule all:
    input:
        ### check if genomic fasta index file exists
        ### needed for rsem and cufflinks
        config["fasta"] + ".fai",
        ### check if fastq files are present otherwise download sra
        ancient(expand(FASTQ_DIR + "{file}", file=Metadata.fq)),
        ### subsample
        get_subsample_pe_output_files,
        get_subsample_se_output_files,
        get_subsample_log_files,
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
        get_multiqc_output_files,
        ### result archive compressed output
        get_result_archive_output_files


#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Execute after the workflow =====                                      │
#└─────────────────────────────────────────────────────────────────────────────┘
onsuccess:
    from datetime import datetime
    print_header(str(datetime.today().strftime("%Y-%m-%d %H:%M:%S")) + 
    " - RNA-Seq workflow finished successfully!")

onerror:
    from datetime import datetime
    print_header(str(datetime.today().strftime("%Y-%m-%d %H:%M:%S")) + 
    " - There was and ERROR running RNA-Seq workflow!")