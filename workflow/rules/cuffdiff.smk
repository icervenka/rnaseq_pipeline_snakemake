OUTDIR = DIFFEXP_OUTDIR + DIFFEXP_ANALYSIS
labels_str, files_str = get_cuffdiff_data()

# TODO add cuffnorm
# TODO switch to conda environment
rule diffexp:
    input:
        gtf = rules.merge.output,
        fasta = config['fasta']
    output:
        OUTDIR + "gene_exp.diff"
    threads:
        config['threads']
    run:
        shell(
            "export PS1=; "
            "source /usr/local/bin/miniconda3/etc/profile.d/conda.sh; "
            "conda activate tophat; "
            "cuffdiff "
            "-p {threads} "
            "-o {OUTDIR} "
            "-b {input.fasta} "
            "-u {input.gtf} "
            "--labels {labels_str} "
            "{files_str} "
        )

rule all_diffexp:
    input:
        rules.diffexp.output
    output:
        touch(LOG_DIR + DIFFEXP_ANALYSIS + "diffexp.completed")
