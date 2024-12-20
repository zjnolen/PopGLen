# Pairwise individual relatedness with NGSrelate and IBSrelate


localrules:
    compile_kinship_stats_sfs,
    ngsrelate_freqbased_merge,
    ngsrelate_freqbased_merge,


rule est_kinship_stats_sfs:
    """
    Uses the equations from Waples et al. 2019, MolEcol to estimate R0, R1, and
    KING-robust kinship between all sample pairings.
    """
    input:
        sfs="results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{ind1}-{ind2}{dp}_{sites}-filts.sfs",
    output:
        "results/datasets/{dataset}/analyses/kinship/ibsrelate_sfs/{dataset}.{ref}_{ind1}-{ind2}{dp}_{sites}-filts.kinship",
    log:
        "logs/{dataset}/kinship/ibsrelate_sfs/{dataset}.{ref}_{ind1}-{ind2}{dp}_{sites}-filts_kinship.log",
    benchmark:
        "benchmarks/{dataset}/kinship/ibsrelate_sfs/{dataset}.{ref}_{ind1}-{ind2}{dp}_{sites}-filts_kinship.log"
    wildcard_constraints:
        ind1="|".join([i for i in samples.index.tolist()]),
        ind2="|".join([i for i in samples.index.tolist()]),
    container:
        r_container
    resources:
        runtime="1h",
    group:
        "sfs-ibsrelate"
    script:
        "../scripts/kinship_sfs.R"


rule compile_kinship_stats_sfs:
    """
    Compiles kinship stats for all pairs into a single table.
    """
    input:
        get_kinship,
    output:
        "results/datasets/{dataset}/analyses/kinship/ibsrelate_sfs/{dataset}.{ref}_all{dp}_{sites}-filts.kinship",
    log:
        "logs/{dataset}/kinship/ibsrelate_sfs/{dataset}.{ref}_all{dp}_{sites}-filts_compile-stats.log",
    benchmark:
        "benchmarks/{dataset}/kinship/ibsrelate_sfs/{dataset}.{ref}_all{dp}_{sites}-filts_compile-stats.log"
    container:
        shell_container
    resources:
        runtime="15m",
    group:
        "sfs-ibsrelate"
    shell:
        """
        (printf "ind1\tind2\tR0\tR1\tKING\n" > {output}
        cat {input} >> {output}) 2> {log}
        """


rule doGlf1_ibsrelate:
    """
    Generate binary format genotype likelihoods for all filtered positions for
    all samples to run IBSrelate IBS-based version
    """
    input:
        unpack(filt_depth),
        bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
    output:
        glf=temp(
            "results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.glf.gz"
        ),
        arg="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.arg",
    log:
        "logs/{dataset}/angsd/doGlf1/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doGlf1/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    container:
        angsd_container
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: lambda wildcards, attempt: attempt * 2
    shell:
        """
        angsd -doGlf 1 -bam {input.bam} -GL {params.gl_model} -ref {input.ref} \
            -nThreads {threads} {params.extra} -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -sites {input.sites} -rf {input.regions} \
            -out {params.out} &> {log}
        """


rule ibsrelate:
    """
    Run IBS-based version of IBS relate for all pairs of samples
    """
    input:
        "results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.glf.gz",
    output:
        "results/datasets/{dataset}/analyses/kinship/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.ibspair",
    log:
        "logs/{dataset}/angsd/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    container:
        angsd_container
    params:
        nind=get_nind,
        maxsites=config["chunk_size"],
        out=lambda w, output: os.path.splitext(output[0])[0],
    resources:
        runtime="7d",
    threads: lambda wildcards, attempt: attempt * 10
    shell:
        """
        ibs -glf {input} -model 0 -nInd {params.nind} -allpairs 1 \
            -m {params.maxsites} -outFileName {params.out}
        """


rule est_kinship_stats_ibs:
    """
    Uses the equations from Waples et al. 2019, MolEcol to estimate R0, R1, and
    KING-robust kinship between all sample pairings.
    """
    input:
        ibs=expand(
            "results/datasets/{{dataset}}/analyses/kinship/ibsrelate_ibs/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.ibspair",
            chunk=chunklist,
        ),
        inds="results/datasets/{dataset}/poplists/{dataset}_{population}{dp}.indiv.list",
    output:
        "results/datasets/{dataset}/analyses/kinship/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_{sites}-filts.kinship",
    log:
        "logs/{dataset}/kinship/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_{sites}-filts_kinship.log",
    benchmark:
        "benchmarks/{dataset}/kinship/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_{sites}-filts_kinship.log"
    container:
        r_container
    resources:
        runtime="1h",
    group:
        "ibs-ibsrelate"
    script:
        "../scripts/kinship_ibs.R"


rule kinship_table_html:
    """
    Converts kinship table to html for report.
    """
    input:
        "results/datasets/{dataset}/analyses/kinship/ibsrelate_{type}/{dataset}.{ref}_all{dp}_{sites}-filts.kinship",
    output:
        report(
            "results/datasets/{dataset}/analyses/kinship/ibsrelate_{type}/{dataset}.{ref}_all{dp}_{sites}-filts.kinship.html",
            category="02 Relatedness",
            subcategory="IBSrelate - {type}",
            labels=lambda w: {"Filter": "{sites}", **dp_report(w), "Type": "Table"},
        ),
    log:
        "logs/{dataset}/kinship/ibsrelate_{type}/{dataset}.{ref}_all{dp}_{sites}-filts_tsv2html.log",
    benchmark:
        "benchmarks/{dataset}/kinship/ibsrelate_{type}/{dataset}.{ref}_all{dp}_{sites}-filts_tsv2html.log"
    container:
        r_container
    shadow:
        "minimal"
    resources:
        runtime="15m",
    group:
        "ibs-ibsrelate"
    script:
        "../scripts/tsv2html.R"


rule ngsrelate_ibsrelate_only:
    """
    Estimates inbreeding and relatedness measures using NGSrelate. Only the
    coefficients for the IBS relate method will be estimated and kept, which
    does not require allele frequencies.
    """
    input:
        beagle="results/datasets/{dataset}/beagles/{dataset}.{ref}_{population}{dp}_{sites}-filts.beagle.gz",
        inds="results/datasets/{dataset}/poplists/{dataset}_{population}{dp}.indiv.list",
    output:
        ngsrelate=temp(
            "results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_{population}{dp}_{sites}-filts_ngsrelate-nofreq.tsv"
        ),
        ibsrelate="results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_{population}{dp}_{sites}-filts_ibsrelate-nofreq.tsv",
        samples="results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_{population}{dp}_{sites}-filts_samples.ibsrelate-nofreq.list",
    wildcard_constraints:
        population="all",
    log:
        "logs/{dataset}/kinship/ngsrelate/{dataset}.{ref}_{population}{dp}_{sites}-filts.ibsrelate-nofreq.log",
    container:
        ngsrelate_container
    threads: lambda wildcards, attempt: attempt * 4
    params:
        nind=get_nind,
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
    shell:
        r"""
        (nsites=$(zcat {input.beagle} | tail -n +2 | wc -l)
        echo "nsites nind"
        echo $nsites {params.nind}
        cut -f1 {input.inds} | tail -n +2 > {output.samples}
        ngsRelate -G {input.beagle} -n {params.nind} -L $nsites \
            -O {output.ngsrelate} -z {output.samples}
        cut -f3-5,30-35 {output.ngsrelate} > {output.ibsrelate}) &> {log}
        """


rule ngsrelate_freqbased:
    """
    Estimates inbreeding and relatedness measures using NGSrelate. This will use
    allele frequencies, enabling all coefficients, including IBS relate ones,
    all will be kept.
    """
    input:
        beagle="results/datasets/{dataset}/beagles/{dataset}.{ref}_{population}{dp}_{sites}-filts.beagle.gz",
        mafs="results/datasets/{dataset}/beagles/{dataset}.{ref}_{population}{dp}_{sites}-filts.mafs.gz",
        inds="results/datasets/{dataset}/poplists/{dataset}_{population}{dp}.indiv.list",
    output:
        relate="results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_{population}{dp}_{sites}-filts_ngsrelate-freq.tsv",
        freq="results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_{population}{dp}_{sites}-filts_ngsrelate-freq.freqs",
        samples="results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_{population}{dp}_{sites}-filts_samples.list",
    wildcard_constraints:
        population="|".join(pop_list),
    log:
        "logs/{dataset}/kinship/ngsrelate/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    container:
        ngsrelate_container
    threads: lambda wildcards, attempt: attempt * 4
    params:
        nind=get_nind,
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
    shell:
        r"""
        (nsites=$(zcat {input.beagle} | tail -n +2 | wc -l)
        echo "nsites nind"
        echo $nsites {params.nind}
        cut -f1 {input.inds} | tail -n +2 > {output.samples}
        zcat {input.mafs} | cut -f7 | sed 1d > {output.freq}
        ngsRelate -G {input.beagle} -n {params.nind} -L $nsites -O {output.relate} \
            -z {output.samples}) &> {log}
        """


rule ngsrelate_freqbased_merge:
    """
    Merge NGSrelate (freq-based) relatedness estimates for all populations into
    one table for the dataset.
    """
    input:
        relates=expand(
            "results/datasets/{{dataset}}/analyses/kinship/ngsrelate/{{dataset}}.{{ref}}_{population}{{dp}}_{{sites}}-filts_ngsrelate-freq.tsv",
            population=pop_list,
        ),
    output:
        "results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts_ngsrelate-freq.tsv",
    log:
        "logs/{dataset}/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts.ngsrelate-freq.log",
    container:
        shell_container
    resources:
        runtime="15m",
    shell:
        """
        (head -n 1 {input.relates[0]} > {output}
        for i in {input.relates}; do
            tail -n+2 $i >> {output}
        done) 2> {log}
        """


rule ngsrelate_summary:
    """
    Converts NGSrelate table to html.
    """
    input:
        "results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts_{method}.tsv",
    output:
        report(
            "results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts_{method}.html",
            category="02 Relatedness",
            subcategory="NgsRelate",
            labels=lambda w: {
                "Filter": "{sites}",
                **dp_report(w),
                "Method": "{method}",
                "Type": "Table",
            },
        ),
    log:
        "logs/{dataset}/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts_{method}.tsv2html.log",
    container:
        r_container
    shadow:
        "minimal"
    resources:
        runtime="15m",
    script:
        "../scripts/tsv2html.R"
