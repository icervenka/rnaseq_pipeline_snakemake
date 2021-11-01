OUTDIR = DIFFEXP_OUTDIR + DIFFEXP_ANALYSIS
DEGFILES_OUTDIR = OUTDIR + "degfiles/"
REPORTS_OUTDIR = OUTDIR + "reports/"
DDS_OUTDIR = OUTDIR + "dds/"

def get_diffexp_output_files(wildcards):
    return [DDS_OUTDIR + "dds.rds",
        DDS_OUTDIR + "result_array.rds",
        DDS_OUTDIR + "result_array_ids.rds",
        OUTDIR + "analysis_params/config.yaml",
        OUTDIR + "analysis_params/metadata.tsv",
        DEGFILES_OUTDIR + "sample_expression.csv",
        REPORTS_OUTDIR + "report.html",
        REPORTS_OUTDIR + "pca.html",
        REPORTS_OUTDIR + "mds_plot.html"]

rule diffexp_init:
    input:
        count_data="counts/counts.tsv",
        col_data="metadata.tsv"
    output:
        DDS_OUTDIR + "dds.rds"
    params:
        design=config['diffexp']["design"],
        ref_levels=config['diffexp']["reference_levels"],
        min_count=config['diffexp']["gene_min_readcount"],
        contrast_type=config['diffexp']["contrast_type"],
        lfc_shrink=config_extra['diffexp']["deseq_lfc_shrink"],
    script:
        "../scripts/deseq_init.R"

# TODO save file with all the contrasts to read by save and gsea
rule diffexp_results:
    input:
        dds=rules.diffexp_init.output,
        col_data=rules.diffexp_init.input.col_data
    output:
        DDS_OUTDIR + "result_array.rds",
    params:
        contrast_type=config['diffexp']["contrast_type"],
        contrasts=config['diffexp']["contrasts"],
        lfc_shrink=config_extra['diffexp']["deseq_lfc_shrink"],
        fdr=config['diffexp']["fdr"]
    script:
        "../scripts/deseq_results.R"

# solve the issue with diffexp csv files
# currently they are created as s side effect
rule diffexp_save:
    input:
        dds=rules.diffexp_init.output,
        result_array=rules.diffexp_results.output[0]
    output:
        sample_expression=DEGFILES_OUTDIR + "sample_expression.csv",
        result_array_ids=DDS_OUTDIR + "result_array_ids.rds",
    params:
        outdir=DEGFILES_OUTDIR,
        species=config["species"],
        gene_ids_in=config['diffexp']["input_gene_ids"],
        gene_ids_out=config['diffexp']["output_gene_ids"]
    script:
        "../scripts/deseq_save.R"

rule diffexp_report:
    input:
        rules.diffexp_init.output,
        rules.diffexp_results.output[0],
        rules.diffexp_save.output.result_array_ids
    output:
        REPORTS_OUTDIR + "report.html",
        REPORTS_OUTDIR + "mds_plot.html"
    params:
        experiment_name=config["experiment_name"],
        outdir=REPORTS_OUTDIR,
        fdr=config['diffexp']["fdr"],
        report_layout=config["report"]["layout"],
        individual=config["report"]["individual"],
        pca_groups=config["report"]["pca"]["groups"],
        heatmap_mode=config["report"]["sample_heatmap"]["mode"],
        heatmap_n_genes=config["report"]["sample_heatmap"]["n_genes"],
        heatmap_cluster_metric_col=config["report"]["sample_heatmap"]["cluster_metric_col"],
        heatmap_cluster_metric_row=config["report"]["sample_heatmap"]["cluster_metric_row"],
        annotation_heatmap=config["report"]["sample_heatmap"]["annotation"],
        annotation_sts=config["report"]["sample_distances_heatmap"]["annotation"],
        upset_maxgroups=config["report"]["upset"]["max_groups"],
        upset_min_group_size=config["report"]["upset"]["min_group_size"],
        species=config["species"],
        gene_ids_in=config['diffexp']["input_gene_ids"],
        mdplot_group=config["report"]["mdplot_group"]
    script:
        "../scripts/deseq_report.R"

rule diffexp_report_pca:
    input:
        dds=rules.diffexp_init.output[0]
    output:
        REPORTS_OUTDIR + "pca.html"
    params:
        species=config['species'],
        id_type=config['diffexp']["input_gene_ids"],
        group=config['report']["pca"]['groups'],
        ngenes=config['report']["pca"]['ngenes'],
        top_pcas=config['report']["pca"]['top_pcas'],
        top_gene_loadings=config['report']["pca"]['top_gene_loadings'],
        pca2go_ngenes=config['report']["pca"]['pca2go_ngenes'],
        pca2go_loadings_ngenes=config['report']["pca"]['pca2go_loadings_ngenes'],
    script:
        "../scripts/deseq_report_pca.R"

rule diffexp_copy_config:
    input:
        config="config.yaml",
        metadata="metadata.tsv"
    output:
        config=OUTDIR + "analysis_params/config.yaml",
        metadata=OUTDIR + "analysis_params/metadata.tsv"
    shell:
        "cp {input.config} {output.config}; cp {input.metadata} {output.metadata}"
