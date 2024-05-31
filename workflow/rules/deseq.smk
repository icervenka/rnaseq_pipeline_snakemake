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
        dds=RDS_OUTDIR + "dds.rds"
    params:
        design=config['diffexp']["design"],
        ref_levels=config['diffexp']["reference_levels"],
        min_count=config['diffexp']["gene_min_readcount"],
        contrast_type=config['diffexp']["contrast_type"],
        lfc_shrink=config_extra['diffexp']["deseq_lfc_shrink"]
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/deseq_init.R"


# TODO save file with all the contrasts to read by save and gsea
rule diffexp_results:
    input:
        dds=rules.diffexp_init.output.dds,
        col_data=ancient(config["metadata"])
    output:
        result_array=RDS_OUTDIR + "result_array.rds",
    params:
        contrast_type=config['diffexp']["contrast_type"],
        contrasts=config['diffexp']["contrasts"],
        lfc_shrink=config_extra['diffexp']["deseq_lfc_shrink"],
        fdr=config['diffexp']["fdr"]
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/deseq_results.R"


# TODO solve the issue with diffexp csv files
# currently they are created as s side effect
rule diffexp_save:
    input:
        dds=rules.diffexp_init.output.dds,
        result_array=rules.diffexp_results.output.result_array
    output:
        sample_expression=DEGFILES_OUTDIR + "sample_expression.txt",
        result_array_ids=RDS_OUTDIR + "result_array_ids.rds"
    params:
        outdir=DEGFILES_OUTDIR,
        species=config["species"],
        ids_in=config['diffexp']["input_gene_ids"]
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
        report=REPORTS_OUTDIR + "report.html"
    params:
        experiment_name=config["experiment_name"],
        outdir=REPORTS_OUTDIR,
        dir_structure=_ds,
        species=config["species"],
        ids_in=config['diffexp']["input_gene_ids"],
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
        upset_min_group_size=config["report"]["upset"]["min_group_size"]
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/deseq_report.R"


rule diffexp_report_pca:
    input:
        dds=rules.diffexp_init.output.dds
    output:
        report_pca=REPORTS_OUTDIR + "pca.html"
    params:
        outdir=REPORTS_OUTDIR,
        dir_structure=_ds,
        species=config['species'],
        ids_in=config['diffexp']["input_gene_ids"],
        group=config['report']["pca"]['groups'],
        top_pcas=config['report']["pca"]['top_pcas'],
        ngenes=config['report']["pca"]['ngenes'],
        top_gene_loadings=config['report']["pca"]['top_gene_loadings'],
        pca2go_ngenes=config['report']["pca"]['pca2go_ngenes'],
        pca2go_loadings_ngenes=config['report']["pca"]['pca2go_loadings_ngenes']
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/deseq_report_pca.R"


# TODO the expression plots are created based on contrasts and are created as a
## side effect
rule diffexp_report_glimma:
    input:
        dds=rules.diffexp_init.output.dds,
        result_array=rules.diffexp_results.output.result_array
    output:
        report_glimma=REPORTS_OUTDIR + "mds-plot.html"
    params:
        outdir=REPORTS_OUTDIR,
        species=config['species'],
        ids_in=config['diffexp']["input_gene_ids"],
        group=config["report"]["mdplot_group"],
        fdr=config['diffexp']["fdr"]
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/glimma_report.R"


rule diffexp_copy_config:
    input:
        config="config.yaml",
        metadata=ancient(config["metadata"])
    output:
        config=OUTDIR + "analysis_params/config.yaml",
        metadata=OUTDIR + "analysis_params/" + config["metadata"]
    shell:
        "cp {input.config} {output.config}; cp {input.metadata} {output.metadata}"
