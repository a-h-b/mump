include:
    "../Snps/prep.smk"
include:
    "../Snps/mutect.smk"
include:
    "../Snps/filter.smk"

# ctrl per taxonomy

rule ctrl_per_tax:
    input:
        checkBinPerTax
    output:
        touch('Snps/{species}/mutect2.perSample.done')

# master command
rule all_snps:
    input:
        expand("Snps/{species}/mutect2.perSample.done",species=config["snps"]["taxonomy"].replace(" ","_").split(','))
    output:
        touch('status/snps.done')


