#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Output functions =====                                                │
#└─────────────────────────────────────────────────────────────────────────────┘
def get_diffexp_output_files(wildcards):
    return [rules.diffexp_init.output.dds,
        rules.diffexp_results.output.result_array,
        rules.diffexp_save.output.result_array_ids,
        rules.diffexp_save.output.sample_expression,
        rules.diffexp_report.output.report,
        rules.diffexp_report_pca.output.report_pca,
        rules.diffexp_report_glimma.output.report_glimma,
        rules.diffexp_copy_config.output.config,
        rules.diffexp_copy_config.output.metadata]


def get_diffexp_log_files(wildcards):
    return []


#┌─────────────────────────────────────────────────────────────────────────────┐
#│ ===== Rules =====                                                           │
#└─────────────────────────────────────────────────────────────────────────────┘
rule diffexp_init:
    input:
        count_data=rules.counts_to_matrix.output,
        col_data=ancient(config["metadata"])
    output:
        dds=opj(RDS_OUTDIR, "dds.rds")
    params:
        diffexp=config["diffexp"],
        diffexp_extra=config_extra["diffexp"]
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/deseq_init.R"


rule diffexp_results:
    input:
        dds=rules.diffexp_init.output.dds,
        col_data=ancient(config["metadata"])
    output:
        result_array=opj(RDS_OUTDIR, "result_array.rds"),
        contrasts=opj(DEGFILES_OUTDIR, "contrasts.txt"),
    params:
        diffexp=config["diffexp"],
        diffexp_extra=config_extra["diffexp"]
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/deseq_results.R"


rule diffexp_save:
    input:
        dds=rules.diffexp_init.output.dds,
        result_array=rules.diffexp_results.output.result_array
    output:
        sample_expression=opj(DEGFILES_OUTDIR, "sample_expression.txt"),
        result_array_ids=opj(RDS_OUTDIR + "result_array_ids.rds")
    params:
        outdir=DEGFILES_OUTDIR,
        species=config["species"],
        diffexp=config["diffexp"],
        diffexp_extra=config_extra["diffexp"]
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/deseq_save.R"


rule diffexp_report:
    input:
        dds=rules.diffexp_init.output.dds,
        result_array=rules.diffexp_results.output.result_array,
        result_array_ids=rules.diffexp_save.output.result_array_ids
    output:
        report=opj(REPORTS_OUTDIR, "report.html")
    params:
        experiment_name=config["experiment_name"],
        outdir=REPORTS_OUTDIR,
        dir_structure=_ds,
        species=config["species"],
        diffexp=config["diffexp"],
        report=config["report"],
        diffexp_extra=config_extra["diffexp"],
        report_extra=config_extra["report"]
        
        # ids_in=config['diffexp']["input_gene_ids"],
        # fdr=config['diffexp']["fdr"],
        # report_layout=config["report"]["layout"],
        # individual=config["report"]["individual"],
        # pca_groups=config["report"]["pca"]["groups"],
        # heatmap_mode=config["report"]["sample_heatmap"]["mode"],
        # heatmap_n_genes=config["report"]["sample_heatmap"]["n_genes"],
        # heatmap_cluster_metric_col=config["report"]["sample_heatmap"]["cluster_metric_col"],
        # heatmap_cluster_metric_row=config["report"]["sample_heatmap"]["cluster_metric_row"],
        # annotation_heatmap=config["report"]["sample_heatmap"]["annotation"],
        # annotation_sts=config["report"]["sample_distances_heatmap"]["annotation"],
        # upset_maxgroups=config["report"]["upset"]["max_groups"],
        # upset_min_group_size=config["report"]["upset"]["min_group_size"]
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/deseq_report.R"


rule diffexp_report_pca:
    input:
        dds=rules.diffexp_init.output.dds
    output:
        report_pca=opj(REPORTS_OUTDIR, "pca.html")
    params:
        outdir=REPORTS_OUTDIR,
        dir_structure=_ds,
        species=config['species'],
        diffexp=config["diffexp"],
        report=config["report"],
        diffexp_extra=config_extra["diffexp"],
        report_extra=config_extra["report"]

        # ids_in=config['diffexp']["input_gene_ids"],
        # group=config['report']["pca"]['groups'],
        # top_pcas=config['report']["pca"]['top_pcas'],
        # ngenes=config['report']["pca"]['ngenes'],
        # top_gene_loadings=config['report']["pca"]['top_gene_loadings'],
        # pca2go_ngenes=config['report']["pca"]['pca2go_ngenes'],
        # pca2go_loadings_ngenes=config['report']["pca"]['pca2go_loadings_ngenes']
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/deseq_report_pca.R"


# The expression plots are created based on contrasts and are created as a
## side effect
rule diffexp_report_glimma:
    input:
        dds=rules.diffexp_init.output.dds,
        result_array=rules.diffexp_results.output.result_array
    output:
        report_glimma=opj(REPORTS_OUTDIR, "mds-plot.html")
    params:
        outdir=REPORTS_OUTDIR,
        species=config['species'],
        diffexp=config["diffexp"],
        report=config["report"],
        diffexp_extra=config_extra["diffexp"],
        report_extra=config_extra["report"]

        # ids_in=config['diffexp']["input_gene_ids"],
        # fdr=config['diffexp']["fdr"],
        # group=config["report"]["mdplot_group"]
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/glimma_report.R"


rule diffexp_copy_config:
    input:
        config="config.yaml",
        metadata=ancient(config["metadata"])
    output:
        config=opj(ANALYSIS_PARAM_OUTDIR, "config.yaml"),
        metadata=opj(ANALYSIS_PARAM_OUTDIR, config["metadata"])
    shell:
        "cp {input.config} {output.config}; cp {input.metadata} {output.metadata}"
