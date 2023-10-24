

rule mutect2:
    input:
        ref="Snps/{taxon}/{bin}.fa",
        bam="Snps/{taxon}/{bin}.bam",
        bed=getBedForSnps, #"Binning/MetaWrap/bins/{clusterID}/contigs.bed",
        fai="Snps/{taxon}/{bin}.fa.fai",
        dict="Snps/{taxon}/{bin}.dict"
    output:
        "Snps/{taxon}/{bin}.mutect2.vcf.gz"
    params:
        af=config['snps']['af_cutoff']
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    threads:
        getThreads(12)
    log:
        "logs/snps_{taxon}_{bin}_mutect.log"
    message:
        "snps_mutect2: Running mutect2 on {wildcards.taxon} {wildcards.bin}."
    conda:
        ENVDIR + "/mump_snps.yaml"
    shell:
        """
        gatk Mutect2 --af-of-alleles-not-in-resource {params.af} \
         -R {input.ref} -I {input.bam} \
         -L {input.bed} -O {output} &> {log}
        """



