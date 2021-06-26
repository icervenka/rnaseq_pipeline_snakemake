rule kallisto_quant:
    input:
        sample=get_fq,
        index=config['index'],
        gtf=config['gtf']
    output:
        h5=ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".h5",
        tsv=ALIGN_OUTDIR + "{sample}/" + KALLISTO_BAM_NAME + ".tsv"
    params:
        length_mean=extra_config['align']['kallisto_extra']['read_length_mean'],
        length_variance=extra_config['align']['kallisto_extra']['read_length_variance'],
        extra=kallisto_params
    threads:
        config['threads']
    run:
        # TODO: redirect to log
        meta = Metadata
        meta['fq_full'] = FASTQ_DIR + meta['fq']
        fq_meta = meta[meta['fq_full'].isin(input.sample)]
        print(fq_meta)
        # TODO: test if ordering is ok on multilane fastqfiles
        input_arr = fq_meta.sort_values(by=['sample']).groupby(['read'])['fq_full'].apply(lambda x: x.to_list())
        input_arr = [' '.join(x) for x in input_arr.to_list()]
        input_str = " ".join(input_arr)
        print(input_str)

        if(fq_meta['read'].unique().size == 1):
            mode = "--single "
            fragment_info = "-l " + str(params.length_mean) + " -s " + str(params.length_variance)
        else:
            mode = " "
            fragment_info = " "

        stranded = kallisto_stranded(fq_meta['stranded'].iloc[0])
        print(stranded)
        shell(
            "kallisto quant "
            "-b 100 "
            "-t {threads} "
            "--genomebam "
            "-g {input.gtf} "
            "-i {input.index} "
            "-o {ALIGN_OUTDIR}{wildcards.sample} "
            "{mode} "
            "{fragment_info} "
            "{stranded} "
            "{params.extra} "
            "{input_str} "
        )

rule all_align:
    input:
        expand(ALIGN_OUTDIR+"{sample}/"+KALLISTO_BAM_NAME+".h5", sample = Samples),
        expand(ALIGN_OUTDIR+"{sample}/"+KALLISTO_BAM_NAME+".tsv", sample = Samples)
    output:
        touch(LOG_DIR + "align.completed"),
        touch(LOG_DIR + "count.completed")
