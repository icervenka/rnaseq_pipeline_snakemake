rule all_diffexp:
    output:
        touch(LOG_DIR + DIFFEXP_ANALYSIS + "diffexp.completed")
