localrules: prep_dRep_single


checkpoint prep_dRep_single:
    input:
        table="sample_table.tsv",
        bins= expand("{raw_directory}/Binning/selected_DASTool_summary.txt",
                     raw_directory=samples.loc[samples.bins == 1,"path"]),
        indir= expand("{raw_directory}/Binning/selected_DASTool_bins",
                     raw_directory=samples.loc[samples.bins == 1,"path"])
    output:
        bins= directory("Dereplication/dRepSingle_genomes")
    resources:
        runtime = "1:00:00",
        mem = MEMCORE
    threads: 1
    message: "prep_dRep_single: Copying DAStool bins to common directory."
    run:
        shell("mkdir {output}")
        for SAMPLE in samples.loc[samples.bins == 1,"sample"]:
            SAMPLE_DIR = samples.loc[samples['sample'] == SAMPLE,"path"].values[0]
            SAMPLE_TAB = os.path.join(SAMPLE_DIR,"Binning/selected_DASTool_summary.txt")
            selected_bins=pd.read_csv(SAMPLE_TAB,sep="\t")
            for BIN in selected_bins['bin']:
                shell("cp {SAMPLE_DIR}/Binning/selected_DASTool_bins/{BIN}.contigs.fa {output.bins}/{SAMPLE}_bin.{BIN}.contigs.fa")

rule prep_checkM_proteins:
    input:
        indir= lambda wildcards: expand("{raw_directory}/Analysis/annotation",
                raw_directory=samples.path[wildcards.sample])
    output:
        "Dereplication/dRepSingle_proteins_perSample/{sample}.prokka.faa.index"
    params:
        idxdir= "Dereplication/dRepSingle_proteins_perSample"
    resources:
        runtime = "1:00:00",
        mem = MEMCORE
    threads: 1
    message: "prep_checkM_proteins: Copying proteins of {wildcards.sample} to common directory and indexing."
    conda: os.path.join(ENVDIR, "/IMP_bioperl.yaml")
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        mkdir -p {params.idxdir}
        cp {input.indir}/prokka.faa {params.idxdir}/{wildcards.sample}.prokka.faa
        {SRCDIR}/fastaIndex.pl {params.idxdir}/{wildcards.sample}.prokka.faa
        """

#instead of dRep_simple:
# extract ORFs per bin to avoid re-calling genes in checkM
rule prep_checkM_ori_faa:
    input:
        aa_index= "Dereplication/dRepSingle_proteins_perSample/{sample}.prokka.faa.index",
        contigs= "Dereplication/dRepSingle_genomes/{sample}_bin.{bin}.contigs.fa",
        gff= lambda wildcards: expand("{raw_directory}/Analysis/annotation/annotation_CDS_RNA_hmms.gff",
                raw_directory=samples.path[wildcards.sample])
    output:
        genes= "Dereplication/dRepSingle_proteins_perBin/{sample}_bin.{bin}.proteins.faa"
    params:
        aa= lambda wildcards: "Dereplication/dRepSingle_proteins_perSample/" + wildcards.sample + ".prokka.faa"
    resources:
        runtime = "1:00:00",
        mem = MEMCORE
    threads: 1
    message: "prep_checkM_ori_faa: Preparing proteins of {wildcards.sample} bin {wildcards.bin} for checkM."
    conda: os.path.join(ENVDIR, "IMP_bioperl.yaml")
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        {SRCDIR}/fastaExtract_bin.pl {params.aa} {input.contigs} {input.gff} {output.genes}
        """

def aggregate_bin_proteins(wildcards):
    checkpoint_output = checkpoints.prep_dRep_single.get(**wildcards).output[0]
    return expand("Dereplication/dRepSingle_proteins_perBin/{sample_bin}.proteins.faa",
           sample_bin=glob_wildcards(os.path.join(checkpoint_output, "{sample_bin}.contigs.fa")).sample_bin)


# run checkM on all bins 
if config['dereplication']['dRep']['checkM_method'] == "taxonomy_wf":
    CHECKM_SHELL="""
        checkm taxonomy_wf -g --pplacer_threads 1 -t {threads} -x faa --tab_table -f {output} domain Bacteria {params.indir} {params.outdir}
        sed -i "s#^Bin Id#genome#" {output}
    """
else:
    CHECKM_SHELL="""
        checkm lineage_wf -g --pplacer_threads 1 -t {threads} -x faa --tab_table -f {output} {params.indir} {params.outdir}
        sed -i "s#^Bin Id#genome#" {output}
        """
rule checkM_ori_faa:
    input:
        genes= aggregate_bin_proteins
    output:
        checkM= "Dereplication/dRepSingle_checkM.tsv"
    params:
        indir= "Dereplication/dRepSingle_proteins_perBin",
        outdir= "Dereplication",
        checkM_method = config['dereplication']['dRep']['checkM_method']
    threads: getThreads(12)
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    conda: os.path.join(ENVDIR, "IMP_multi.yaml")
    shell:
        CHECKM_SHELL



# use checkM output for dRep steps
rule dRep_simple:
    input:
        checkM="Dereplication/dRepSingle_checkM.tsv",
        dir="Dereplication/dRepSingle_genomes"
    output:
        directory("Dereplication/dRepSingle_out/dereplicated_genomes"),
        expand("Dereplication/dRepSingle_out/data_tables/{file}.csv",
                file=["Chdb","Cdb","Sdb","Wdb","Widb","Ndb"])
    threads: getThreads(BIGCORENO)
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    params:
        out = "Dereplication/dRepSingle_out",
        db = config['db_path'],
        length = config['dereplication']['dRep']['length'],
        completeness = config['dereplication']['dRep']['completeness'],
        contamination = config['dereplication']['dRep']['contamination'],
        MASH_sketch = config['dereplication']['dRep']['MASH_sketch'],
        S_algorithm = config['dereplication']['dRep']['S_algorithm'],
        P_ani = config['dereplication']['dRep']['P_ani'],
        S_ani = config['dereplication']['dRep']['S_ani'],
        SkipMash = config['dereplication']['dRep']['SkipMash'],
        SkipSecondary = config['dereplication']['dRep']['SkipSecondary'],
        cov_thresh = config['dereplication']['dRep']['cov_thresh'],
        coverage_method = config['dereplication']['dRep']['coverage_method'],
        clusterAlg = config['dereplication']['dRep']['clusterAlg'],
        checkM_method = config['dereplication']['dRep']['checkM_method'],
        completeness_weight = config['dereplication']['dRep']['completeness_weight'],
        contamination_weight = config['dereplication']['dRep']['contamination_weight'],
        strain_heterogeneity_weight = config['dereplication']['dRep']['strain_heterogeneity_weight'],
        N50_weight = config['dereplication']['dRep']['N50_weight'],
        size_weight = config['dereplication']['dRep']['size_weight']
    log: "logs/dereplication_dRep_simple.log"
    message: "dRep_simple: Running dRep on bins of input samples."
    conda: os.path.join(ENVDIR, "IMP_multi.yaml")
    shell:
        """
        dRep filter {params.out} -d -p {threads} --length {params.length} --completeness {params.completeness} \
         --contamination {params.contamination} --genomeInfo {input.checkM}  \
         -g {input.dir}/*contigs.fa >> {log} 2>&1
        dRep cluster {params.out} -p {threads} --MASH_sketch {params.MASH_sketch} --S_algorithm {params.S_algorithm} \
         --P_ani {params.P_ani} --S_ani {params.S_ani} {params.SkipMash} \
         {params.SkipSecondary} --cov_thresh {params.cov_thresh} \
         --coverage_method {params.coverage_method} --clusterAlg {params.clusterAlg} >> {log} 2>&1
        dRep choose {params.out} -p {threads} --completeness_weight {params.completeness_weight} \
         --contamination_weight {params.contamination_weight} --strain_heterogeneity_weight {params.strain_heterogeneity_weight} \
         --N50_weight {params.N50_weight} --size_weight {params.size_weight} --checkM_method {params.checkM_method} >> {log} 2>&1
        dRep evaluate {params.out} --evaluate all -p {threads} >> {log} 2>&1
        """
        
