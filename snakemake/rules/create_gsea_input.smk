GSEA_OUT = DIFFEXP_OUTDIR + DIFFEXP_ANALYSIS + GSEA_INPUT_OUTDIR

def get_gsea_output(wildcards):
    return [
        GSEA_OUT + config['experiment_name'] + "gsea_expression.gct",
        GSEA_OUT + config['experiment_name'] + ".cls"
    ] + 
    expand(GSEA_OUT + "{contrast}.rnk", contrast=get_contrasts)

rule expression2gsea:
    input:
        get_sample_expression
    output:
        expression=GSEA_OUT + config['experiment_name'] + "gsea_expression.gct"
        classes=GSEA_OUT + config['experiment_name'] + ".cls"
    threads:
        1
    script:
        "../scripts/expression2gsea.R"

rule expression2gsea_preranked:
    input:
        get_contrast_files
    output:
        GSEA_OUT + "{contrast}.rnk"
    params:
        ranking=config['gsea']['ranking']
    threads:
        1
    script:
        "../scripts/expression2gsea_preranked.R"
