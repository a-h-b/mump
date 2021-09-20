include:
    "../Collate/features.smk"
include:
    "../Collate/taxonomy.smk"
include:
    "../Collate/stats.smk"

# master command
rule all_collate:
    input:
        expand("status/{task}.collated",task=config["collection"]["results"].split())
    output:
        touch('status/collate.done')


