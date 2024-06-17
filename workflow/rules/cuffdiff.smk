def get_diffexp_output_files(wildcards):
    return (
        rules.diffexp.output + 
        rules.scatter_diffexp.output
    )

def get_diffexp_log_files(wildcards):
    return [rules.diffexp.log]


rule make_cuffdiff_gtf:
    input:
        lambda wildcards: get_gtf()
    output:
        gtf=temp(opj(CUFFCOMPARE_OUTDIR, "cuffdiff" + CUFFCOMPARE_GTF_FILE)),
        other=temp(expand(opj(CUFFCOMPARE_OUTDIR, "cuffdiff" + "{files}"), files=CUFFCOMPARE_FILES))
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
        gtf=rules.make_cuffdiff_gtf.output.gtf,
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


rule scatter_diffexp:
    input:
        lambda wildcards: [ x for x in rules.diffexp.output if CUFFDIFF_GENE_DEG_FILE in x ]
    output:
        opj(DEGFILES_OUTDIR, "cuffdiff_contrasts.txt")
    params:
        outdir=DEGFILES_OUTDIR,
        fdr=config["diffexp"]["fdr"]
    threads:
        1
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/cuffdiff_separate.R"


include: "copy_config.smk"