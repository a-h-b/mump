binctrl = []
if config['dereplication']['simple_derep']:
    include:
        "../../Dereplicate/dRep.simple.smk"
    include:
            "../../Dereplicate/dRep_sum.simple.smk"
    binctrl.append("Dereplication/dereplicated_bin_stats.tsv")
if config['dereplication']['cross_mapping_rebinning']['do']:
    include:
        "../../Dereplicate/cross_mapping.smk"
    include:
        "../../Dereplicate/MetaBAT.smk"
    include:
        "../../Dereplicate/MaxBin.smk"
    include:
        "../../Dereplicate/dastool.smk"
    include:
        "../../Dereplicate/grid.smk"
    include:
        "../../Dereplicate/dRep.multi.smk"
    binctrl.append("Dereplication/rebinned_bin_stats.tsv")
    if config['dereplication']['GTDBtk']:
        include:
            "../../Dereplicate/gtdbtk_multi.smk"
    else:
        include:
            "../../Dereplicate/dRep_sum.multi.smk"


# master command
rule all_derep:
    input:
        binctrl
    output:
        touch('status/dereplicate.done')


