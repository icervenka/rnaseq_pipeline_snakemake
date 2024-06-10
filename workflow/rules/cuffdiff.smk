def get_diffexp_output_files(wildcards):
    return (
        rules.diffexp.output
    )

def get_diffexp_log_files(wildcards):
    return [rules.diffexp.log]


# gtf has to have certain attributes
# Note: If an arbitrary GTF/GFF3 file is used as input (instead of the 
# .combined.gtf file produced by Cuffcompare), these attributes will not be 
# present, but Cuffcompare can still be used to obtain these attributes with a command like this:
# cuffcompare -s /path/to/genome_seqs.fa -CG -r annotation.gtf annotation.gtf

from pprint import pprint
pprint(vars(rules))
print(print([rule.name for rule in workflow.rules]))

def has_rule(rule):
    rule_names = [rule.name for rule in workflow.rules]
    if rule in rule_names:
        return True
    else:
        return False
    

# TODO maybe put into function
if has_rule("cuffmerge"):
    def get_gtf_for_cuffdiff(wildcards):
        return rules.cuffmerge.output.merged_gtf
elif has_rule("merge"):
    def get_gtf_for_cuffdiff(wildcards):
        return rules.merge.output.merged_gtf
else:
    def get_gtf_for_cuffdiff(wildcards):
        return config["gtf"]


rule make_cuffdiff_gtf:
    input:
        get_gtf_for_cuffdiff
    output:
        temp(opj(COUNT_OUTDIR, "cuffcompare", "cuffcmp.combined.gtf"))
    params:
        outprefix=opj(COUNT_OUTDIR, "cuffcompare", "cuffcmp"),
        fasta=config["fasta"]
    threads:
        config["threads"]
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    shell:
        """
        cuffcompare -s {params.fasta} -o {params.outprefix} -CG {input}
        """


# function for arranging files for input output
rule diffexp:
    input:
        bam=expand(rules.align_out.output, sample=Samples),
        gtf=rules.make_cuffdiff_gtf.output,
    output:
        expand(opj(DEGFILES_OUTDIR, "{file}"), file=CUFFDIFF_DIFFEXP_FILES)
    params:
        outdir=DEGFILES_OUTDIR,
        metadata=Metadata,
        samples=Samples,
        # pass directory structure and tool filenames into the wrapper
        ds=_ds,
        tf=_tf,
        diffexp=config["diffexp"],
        fasta=config["fasta"]
    log:
        log=expand(opj(DIFFEXP_LOG_OUTDIR, "{log}"), log=CUFFDIFF_LOG_FILES)
    threads:
        config['threads']
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    script:
        "../scripts/cuffdiff_wrapper.py"


    # shell:
    #     """
    #     cuffdiff \
    #     --no-update-check
    #     -q \
    #     -p {threads} \
    #     -o {params.outdir} \
    #     --frag-bias-correct {params.fasta} \
    #     --FDR
    #     -u {input.gtf} \
    #     --labels {params.labels} \
    #     {params.files} \
    #     > {log} 2>&1
    #     """


# rule separate_diffexp:
#     input:
#         ""
#     output:
#         ""
#     params:
#         ""
#     script:
#         ""

