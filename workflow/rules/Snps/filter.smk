rule snps_filter:
    input:
        "Snps/{taxon}/{bin}.mutect2.vcf.gz"
    output:
        "Snps/{taxon}/{bin}.mutect2.filtered.vcf.gz"
    params:
        af=config['snps']['af_cutoff']
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    log: "logs/snps_filter.{taxon}.{bin}.log"
    message: 'snps_filter: filtering mutect output of {wildcards.taxon} {wildcards.bin}'
    threads: 1
    conda: 
        os.path.join(ENVDIR, "mump_snps.yaml")
    script:
        os.path.join(SRCDIR, "snps_filtering_snakemake.R")



