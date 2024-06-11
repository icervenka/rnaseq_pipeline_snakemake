def get_diffexp_output_files(wildcards):
    return (
        rules.diffexp.output
    )

def get_diffexp_log_files(wildcards):
    return [rules.diffexp.log]

 
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
        temp(expand(opj(CUFFCOMPARE_OUTDIR, "cuffdiff" + "{files}"), files=CUFFCOMPARE_NAMES))
    params:
        outprefix=opj(CUFFCOMPARE_OUTDIR, "cuffdiff"),
        fasta=config["fasta"]
    threads:
        config["threads"]
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    shell:
        """
        cuffcompare -s {params.fasta} -o {params.outprefix} -CG {input}
        """


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


# rule separate_diffexp:
#     input:
#         ""
#     output:
#         ""
#     params:
#         ""
#     script:
#         ""

include: "copy_config.smk"