def get_count_output_files(wildcards):
    return (
        expand(rules.cufflinks.output, sample=Samples) +
        rules.cuffmerge.output +
        [COUNT_OUTDIR + "abundances.cxb"]
    )

def get_count_log_files(wildcards):
    return []

rule cufflinks:
    input:
        bam=rules.align_out.output,
        gtf=config["gtf"],
        fasta=config["fasta"]
    output:
        COUNT_OUTDIR + "{sample}/" + CUFFLINKS_GTF_FILE
    params:
        outdir=COUNT_OUTDIR,
        stranded="", #cufflink_stranded,
        standard="", #cufflink_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    threads:
        config["threads"]
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    shell:
        """
        cufflinks \
        --frag-bias-correct {input.fasta} \
        {params.stranded} \
        {params.standard} \
        {params.extra} \
        -o {params.outdir}{wildcards.sample} \
        {input.bam}
        """

rule gather_gtf:
    input:
        expand(rules.cufflinks.output, sample=Samples)
    output:
        temp(COUNT_LOG_OUTDIR + "gtf_assembly.txt")
    shell:
        "printf '%s\n' {input} >> {output}"

# TODO check if gtf should be specified also for denovo assembly
rule cuffmerge:
    input:
        gathered=rules.gather_gtf.output,
        gtf=config['gtf'],
        fasta=config["fasta"]
    output:
        COUNT_OUTDIR + CUFFLINKS_MERGED_FILE
    params:
        outdir=COUNT_OUTDIR,
        stranded="", #cufflink_stranded,
        standard="", #cufflink_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    threads:
        config['threads']
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    shell:
        """
        cuffmerge \
        -p {threads} \
        -g {input.gtf} \
        -s {input.fasta} \
        {params.stranded} \
        {params.standard} \
        {params.extra} \
        -o {params.outdir} \
        {input.gathered}
        """

rule cuffcompare:
    input:
        expand(rules.cufflinks.output, sample=Samples)
    output:
        COUNT_OUTDIR + "cuffcompare/" + "cuffcmp.loci"
    params:
        outdir=COUNT_OUTDIR + "cuffcompare/",
        gtf=config["gtf"]
    threads:
        config['threads']
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    shell:
        "cuffcompare {input} -o {output}"


rule cuffquant:
    input:
        bam=expand(rules.align_out.output, sample=Samples),
        gtf=rules.cuffmerge.output,
        fasta=config["fasta"]
    output:
        COUNT_OUTDIR + "abundances.cxb"
    params:
        outdir=COUNT_OUTDIR,
        stranded="", #cufflink_stranded,
        standard="", #cufflink_params,
        extra=has_extra_config(config["count"]["extra"], config_extra["count"])
    threads:
        config['threads']
    conda:
        CONDA_COUNT_CUFFLINKS_ENV
    shell:
        """
        cuffquant \
        --frag-bias-correct {input.fasta}
        {params.stranded} \
        {params.standard} \
        {params.extra} \
        -p {threads} \
        {input.gtf} \
        {input.bam}
        """
