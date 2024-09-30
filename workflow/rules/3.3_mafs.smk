rule link_mafs:
    """
    The population inferred major allele MAFs are the same as those made during
    population level SNP calling and can just be linked.
    """
    input:
        "results/datasets/{dataset}/beagles/{dataset}.{ref}_{population}{dp}_{sites}-filts.mafs.gz",
    output:
        "results/datasets/{dataset}/mafs/{dataset}.{ref}_{population}{dp}_{sites}-filts.pop-maj.mafs.gz",
    log:
        "logs/datasets/{dataset}/link_mafs/{dataset}.{ref}_{population}{dp}_{sites}-filts.pop-maj.mafs.gz",
    wildcard_constraints:
        population="|".join(pop_list),
    container:
        shell_container
    shell:
        """
        ln {input} {output} 2> {log}
        """


rule angsd_doMaf:
    """
    Dataset inferred major alleles require calculating MAFs using the SNPs
    called across the dataset. Major and minor allele will match in across all
    populations in this file (this will already be true though, for ref or anc
    assigned major alleles). Another important difference from the population
    level files is that it keeps all sites that are variable in the dataset,
    even if they are not within the population, so these files can be used
    for comparing allele frequencies between populations and tools that do that,
    such as treemix.
    """
    input:
        unpack(get_anc_ref),
        sites="results/datasets/{dataset}/filters/snps/{dataset}.{ref}_all{dp}_{sites}-filts_snps.sites",
        idx="results/datasets/{dataset}/filters/snps/{dataset}.{ref}_all{dp}_{sites}-filts_snps.sites",
        bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        reffai="results/ref/{ref}/{ref}.fa.fai",
    output:
        maf="results/datasets/{dataset}/mafs/{dataset}.{ref}_{population}{dp}_{sites}-filts.dataset-maj.mafs.gz",
        arg="results/datasets/{dataset}/mafs/chunks/{dataset}.{ref}_{population}{dp}_{sites}-filts.dataset-maj.arg",
    log:
        "logs/{dataset}/angsd/doMaf/{dataset}.{ref}_{population}{dp}_{sites}-filts.dataset-maj.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doMaf/{dataset}.{ref}_{population}{dp}_{sites}-filts.dataset-maj.log"
    container:
        angsd_container
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        extra_beagle=config["params"]["angsd"]["extra_beagle"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        maf=config["params"]["angsd"]["domaf"],
        counts=get_docounts,
        trans=get_trans,
        nind=get_nind,
        minind=get_minind,
        mininddp=config["params"]["angsd"]["mindepthind"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        angsd -bam {input.bam} -GL {params.gl_model} -ref {input.ref} \
            -doMajorMinor 3 -doMaf {params.maf} -nThreads {threads} \
            {params.extra} -minMapQ {params.mapQ} -minQ {params.baseQ} \
            -sites {input.sites} -anc {input.anc} {params.extra_beagle} \
            {params.minind} -setMinDepthInd {params.mininddp} {params.counts} \
            -rmTrans {params.trans} -out {params.out} &> {log}
        """
