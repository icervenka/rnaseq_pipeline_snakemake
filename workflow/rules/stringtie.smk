if config["count"]["transcript_assembly"] in ["guided", "guided", "denovo", "Denovo"]:
    include: "stringtie_assembly.smk"
else:
    include: "stringtie_abundance.smk"
