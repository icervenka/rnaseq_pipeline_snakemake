def get_align_output_files(wildcards):
    return
        expand(ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".h5", sample = Samples) +
        expand(ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".tsv", sample = Samples) +
        expand(ALIGN_OUTDIR + "{sample}/" + "pseudoalignments.bam", sample = Samples) +
        expand(ALIGN_OUTDIR + "{sample}/" + "pseudoalignments.bam.bai", sample = Samples)

def get_align_log_files(wildcards):
    return expand(ALIGN_LOG_OUTDIR + "{sample}/" + "run_info.json", sample = Samples)

rule kallisto_quant:
    input:
        sample=get_fq,
        index=config['index'],
        gtf=config['gtf']
    output:
        h5=ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".h5",
        tsv=ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".tsv",
        log=ALIGN_OUTDIR + "{sample}/" + "run_info.json"
    params:
        stranded=kallisto_stranded,
        length_mean=config_extra['align']['kallisto_extra']['read_length_mean'],
        length_variance=config_extra['align']['kallisto_extra']['read_length_variance'],
        extra=config_extra['align']['kallisto_extra']['extra']
    threads:
        config['threads']
    run:
        # # TODO: redirect to log
        # meta = Metadata
        # meta['fq_full'] = FASTQ_DIR + meta['fq']
        # fq_meta = meta[meta['fq_full'].isin(input.sample)]
        # print(fq_meta)
        # # TODO: test if ordering is ok on multilane fastqfiles
        # input_arr = fq_meta.sort_values(by=['sample']).groupby(['read'])['fq_full'].apply(lambda x: x.to_list())
        # input_arr = [' '.join(x) for x in input_arr.to_list()]
        # input_str = " ".join(input_arr)
        # print(input_str)

        input_arr = arrange_fq_for_align(input.sample, Metadata, FASTQ_DIR)

        if len(input_arr) == 1:
            input_str = "--single " + " ".join(input_arr[0])
            fragment_info = "-l " + str(params.length_mean) + " -s " + str(params.length_variance)
        elif len(input_arr) == 2:
            input_str = " ".join(zip(input_arr[0], input_arr[1]))
            fragment_info = ""
        else:
            raise ValueError("Too many read types.")


        # if(fq_meta['read'].unique().size == 1):
        #     mode = "--single "
        #     fragment_info = "-l " + str(params.length_mean) + " -s " + str(params.length_variance)
        # else:
        #     mode = " "
        #     fragment_info = " "

        # stranded = kallisto_stranded(fq_meta['stranded'].iloc[0])
        # print(stranded)
        shell(
            "kallisto quant "
            "-b 100 "
            "-t {threads} "
            "--genomebam "
            "-g {input.gtf} "
            "-i {input.index} "
            "-o {ALIGN_OUTDIR}{wildcards.sample} "
            "{fragment_info} "
            "{params.stranded} "
            "{params.extra} "
            "{input_str} "
        )

rule move_align_log:
    input:
        rules.align.output.log
    output:
        ALIGN_LOG_OUTDIR + "{sample}/" + "run_info.json"
    shell:
        "cp {input} {output}"

# rule all_align:
#     input:
#         expand(ALIGN_OUTDIR+"{sample}/"+KALLISTO_BAM_NAME+".h5", sample = Samples),
#         expand(ALIGN_OUTDIR+"{sample}/"+KALLISTO_BAM_NAME+".tsv", sample = Samples)
#     output:
#         touch(LOG_DIR + "align.completed"),
#         touch(LOG_DIR + "count.completed")
