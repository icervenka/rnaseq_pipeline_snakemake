def get_align_log_files(wildcards):
    return []


def get_bam_index_files(wildcards):
    return []


if "skip_count" == pipeline["diffexp"]:
    def get_align_output_files(wildcards):
        return []
else:
    def get_align_output_files(wildcards):
        return expand(rules.align_out.output, sample=Samples)

    if "sleuth" == pipeline["diffexp"]:
        if "h5" in present_ext:
            COMMON_ALIGN_FILE = KALLISTO_QUANT_H5_FILE
        else:
            COMMON_ALIGN_FILE = SALMON_QUANT_FILE
    else:
        COMMON_ALIGN_FILE = COMMON_BAM_FILE

    rule align_out:
        input:
            glob_wildcards(opj(ALIGN_OUTDIR, "{sample}", "{file}"))
        output:
            expand(opj(ALIGN_OUTDIR, "{{sample}}", "{file}"), file = COMMON_ALIGN_FILE)
        shell:
            "mv {input} {output}"
