rule all_count:
    output:
        touch(LOG_DIR + "count.completed")
