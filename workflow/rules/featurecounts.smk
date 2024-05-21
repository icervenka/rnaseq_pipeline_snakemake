if all_paired or all_single:
    include: "featurecounts_multi.smk"
else:
    include: "featurecounts_single.smk"
