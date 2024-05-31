def get_align_log_files(wildcards):
    return []


def get_bam_index_files(wildcards):
    return []


if "skip_count" in config['pipeline']:
    def get_align_output_files(wildcards):
        return []
else:
    def get_align_output_files(wildcards):
        return expand(rules.align_out.output, sample=Samples)

    if "sleuth" in config['pipeline']:
        if "h5" in present_ext:
            COMMON_ALIGN_NAME = KALLISTO_QUANT_NAME
            ext = "h5"
        else:
            ext = "sf"
            COMMON_ALIGN_NAME = SALMON_QUANT_NAME
    else:
        COMMON_ALIGN_NAME = COMMON_BAM_NAME
        ext = "bam"

    rule align_out:
        input:
            glob_wildcards(opj(ALIGN_OUTDIR, "{sample}", "{file}" + "." + ext))
        output:
            expand(opj(ALIGN_OUTDIR, "{{sample}}", "{file}" + "." + ext), file = COMMON_ALIGN_NAME)
        shell:
            "mv {input} {output}"
