#  def get_stringtie_processed_output_files(wildcards):
#     return (
#         [rules.counts_to_matrix.output.counts_gene,
#         rules.counts_to_matrix.output.counts_transcript, 
#         rules.gather_data.output.tpm,
#         rules.gather_data.output.fpkm,
#         rules.gather_data.output.samples_combined] 
#     )

def get_stringtie_processed_output_files(wildcards):
    return [rules.gather_data.output.tpm,
        rules.gather_data.output.fpkm,
        rules.gather_data.output.samples_combined,
        rules.counts_to_matrix.output[0]]

#  TODO https://github.com/maplexuci/stringtie_gene_id_replacement
# https://www.reddit.com/r/bioinformatics/comments/6xd4py/anyone_have_a_good_script_for_turning_mstrg_gene/

rule counts_to_matrix:
    input:
        files=get_tximport_files,
        gtf=lambda wildcards: get_gtf()
    output:
        opj(COUNT_OUTDIR, COMMON_COUNT_NAME)
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


# TODO remove the ensembl gene version from count matrix
# TODO replace with tximport
# TODO comprare with tximport
##  remove last row with colsums
# rule counts_to_matrix:
#     input:
#         expand(rules.count.output.counts, sample=Samples)
#     output:
#         counts_gene=opj(COUNT_OUTDIR, COMMON_COUNT_NAME,)
#         counts_transcript=opj(COUNT_OUTDIR, COMMON_TRANSCRIPT_COUNT_NAME)
#     params:
#         input_dir=COUNT_OUTDIR,
#         length=75
#     threads:
#         1
#     conda:
#         CONDA_COUNT_GENERAL_ENV
#     shell:
#         """
#         python3 workflow/scripts/prepDE.py \
#         -l {params.length} \
#         -i {params.input_dir} \
#         -g {output.counts_gene} \
#         -t {output.counts_transcript} \
#         """


rule gather_data:
    input:
        expand(rules.count.output.counts, sample=Samples)
    output:
        tpm=opj(COUNT_OUTDIR, STRINGTIE_TPM_FILE),
        fpkm=opj(COUNT_OUTDIR, STRINGTIE_FPKM_FILE),
        samples_combined=opj(COUNT_OUTDIR, STRINGTIE_COMBINED_FILE)
    params:
        samples=expand("{sample}", sample=Samples)
    conda:
        CONDA_R_GENERAL_ENV
    script:
        "../scripts/stringtie_gather.R"