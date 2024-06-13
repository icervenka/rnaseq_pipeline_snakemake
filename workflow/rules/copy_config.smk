rule diffexp_copy_config:
    input:
        config="config.yaml",
        config_extra="config_extra.yaml",
        metadata=ancient(config["metadata"])
    output:
        config=opj(ANALYSIS_PARAM_OUTDIR, "config.yaml"),
        config_extra=opj(ANALYSIS_PARAM_OUTDIR, "config_extra.yaml"),
        metadata=opj(ANALYSIS_PARAM_OUTDIR, config["metadata"])
    shell:
        """
        cp {input.config} {output.config}; 
        cp {input.config_extra} {output.config_extra}; 
        cp {input.metadata} {output.metadata}
        """