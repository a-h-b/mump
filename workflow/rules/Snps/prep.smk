checkpoint prep_bins_snps:
    input:
        table="sample_table.tsv",
        bins= expand("{raw_directory}/Stats/all_bin_stats.tsv",
                     raw_directory=samples.loc[samples.bins == 1,"path"]),
        indir= expand("{raw_directory}/Binning/MetaWrap/bins",
                     raw_directory=samples.loc[samples.bins == 1,"path"])
    output:
        bins= expand("Snps/{taxon}/prepped",taxon=config["snps"]["taxonomy"].replace(" ","_").split(','))
    resources:
        runtime = "1:00:00",
        mem = MEMCORE
    threads: 1
    message: "prep_bins_snps: Copying bins to per species directory."
    run:
        for spec in config["snps"]["taxonomy"].split(','):
            spec_outdir = os.path.join("Snps",spec.replace(" ","_"))
            print(spec)
            os.makedirs(spec_outdir, exist_ok=True)
            for SAMPLE in samples.loc[samples.tax_bins == 1,"sample"]:
                SAMPLE_DIR = samples.loc[samples['sample'] == SAMPLE,"path"].values[0]
                SAMPLE_TAB = os.path.join(SAMPLE_DIR,"Stats/all_bin_stats.tsv")
                SAMPLE_bins=pd.read_csv(SAMPLE_TAB,sep="\t")
                respec = ".*" + spec + "$"
                #print(respec)
                #print(SAMPLE_bins.classification)
                #print(SAMPLE_bins[SAMPLE_bins.classification.str.match(respec)])
                for BIN in SAMPLE_bins.loc[SAMPLE_bins.classification.str.match(respec),"bin"]:
                    shell("cp {SAMPLE_DIR}/Binning/MetaWrap/bins/{BIN}.fasta {spec_outdir}/{SAMPLE}_{BIN}.fa")
                    print(SAMPLE + " : " + BIN)
            shell("touch {spec_outdir}/prepped")


rule contig_fasta_indexing:
    input:
        "Snps/{taxon}/{bin}.fa" 
    output:
        "Snps/{taxon}/{bin}.fa.fai"
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
    conda: ENVDIR + "/mump_snps.yaml"
    log: "logs/snps_{taxon}_{bin}_fasta_indexing.log"
    message: "snps_fasta_indexing: Indexing {wildcards.bin} for {wildcards.taxon}."
    shell:
        """
        samtools faidx {input[0]} > {log} 2>&1
        """

rule contig_dictionary:
    input:
        "Snps/{taxon}/{bin}.fa"    
    output:
        "Snps/{taxon}/{bin}.dict"
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
    conda: ENVDIR + "/mump_snps.yaml"
    log: "logs/snps_{taxon}_{bin}_dictionary.log"
    message: "snps_contig_dictionary: GATK dictionary {wildcards.bin} for {wildcards.taxon}."
    shell:
        """
        gatk CreateSequenceDictionary -R {input} -O {output} > {log} 2>&1
        """


rule prepare_snpBams:
    input:
        ass_bam=getBamForSnps, #"Assembly/mg.reads.sorted.bam",
        bin_bed=getBedForSnps #"Binning/MetaWrap/bins/{clusterID}/contigs.bed"
    output:
        bam="Snps/{taxon}/{bin}.bam",
        bai="Snps/{taxon}/{bin}.bam.bai"
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads:
        1
    log:
        "logs/snps_{taxon}_{bin}_indexbam.log"
    message:
        "snps_prepare_snpBams: Indexing bam file for {wildcards.taxon} {wildcards.bin}."
    conda:
        ENVDIR + "/IMP_mapping.yaml"
    shell:
        """
        samtools view -b -F 3844 -L {input.bin_bed} {input.ass_bam}  2>> {log} | \
          samtools sort -o {output.bam} - > {log} 2>&1
        samtools index {output.bam}
        """
