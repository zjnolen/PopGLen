# Rules for calculating genotype likelihoods from sequencing data using ANGSD


localrules:
    angsd_makeBamlist,
    popfile,


rule angsd_makeBamlist:
    """
    Create list of bam files to be input into calculation.
    """
    input:
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
    output:
        "results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
    log:
        "logs/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.log",
    container:
        shell_container
    params:
        sampord=lambda w, input: f"{input.bams}",
    resources:
        runtime="15m",
    shell:
        """
        (
        echo "BAM order: {params.sampord}" > {log}
        readlink -f {input.bams} > {output} 2>> {log}
        truncate -s -1 {output}) 2>> {log}
        """


rule popfile:
    """
    Create list of samples per population to preserve order of samples as sample names 
    are not propagated into downstream files.
    """
    output:
        inds="results/datasets/{dataset}/poplists/{dataset}_{population}{dp}.indiv.list",
    log:
        "logs/{dataset}/poplists/{dataset}_{population}{dp}_makelist.log",
    container:
        pandas_container
    params:
        samplelist=samples,
        inds=get_popfile_inds,
    resources:
        runtime="15m",
    script:
        "../scripts/make_popfile.py"


rule angsd_sites_index:
    """
    Index sites files used for limiting calculations to specific positions.
    """
    input:
        "results/datasets/{dataset}/filters/{prefix}.sites",
    output:
        multiext("results/datasets/{dataset}/filters/{prefix}.sites", ".idx", ".bin"),
    log:
        "logs/{dataset}/angsd/sites_index/{prefix}.log",
    benchmark:
        "benchmarks/{dataset}/angsd/sites_index/{prefix}.log"
    container:
        angsd_container
    resources:
        runtime="1h",
    shell:
        """
        angsd sites index {input} 2> {log}
        """
