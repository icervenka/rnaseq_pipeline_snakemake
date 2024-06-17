
def get_counts_to_matrix_output_files(wildcards):
    return rules.counts_to_matrix.output


# https://github.com/maplexuci/stringtie_gene_id_replacement
# https://www.reddit.com/r/bioinformatics/comments/6xd4py/anyone_have_a_good_script_for_turning_mstrg_gene/
rule counts_to_matrix:
    input:
        files=get_tximport_files,
        gtf=lambda wildcards: get_gtf() 
    output:
        opj(COUNT_OUTDIR, COMMON_COUNT_FILE)
    params:
        metadata=config["metadata"],
        samples=Samples,
        species=config["species"],
        pipeline=pipeline.values(),
        tool=config_extra["tximport_count_matrix"]["tool"],
        ids_in=config["diffexp"]["input_gene_ids"],
        out_type="counts",
        extra=config_extra["tximport_count_matrix"]["extra"]
    threads:
        1
    conda:
        CONDA_DIFFEXP_GENERAL_ENV
    script:
        "../scripts/tximport.R"