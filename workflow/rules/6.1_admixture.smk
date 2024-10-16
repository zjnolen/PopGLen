# Rules for performing and evaluating individual admixture analyses


rule ngsAdmix:
    """
    Performs individual admixture analyses for a given K, performing several replicates 
    and preserving the replicate with the highest likelihood. After $minreps, 
    convergence is assessed and replicates cease when convergence is reached or $reps
    replicates have been performed.
    """
    input:
        beagle=expand(
            "results/datasets/{{dataset}}/beagles/pruned/{{dataset}}.{{ref}}_{{population}}{{dp}}_{{sites}}-filts.pruned_maxkbdist-{maxkb}_minr2-{r2}.beagle.gz",
            maxkb=config["params"]["ngsld"]["max_kb_dist_pruning_dataset"],
            r2=config["params"]["ngsld"]["pruning_min-weight_dataset"],
        ),
    output:
        qopt="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.qopt",
        fopt="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.fopt.gz",
        log="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}_optimization_wrapper.log",
    log:
        "logs/{dataset}/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.log",
    benchmark:
        "benchmarks/{dataset}/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.log"
    container:
        angsd_container
    params:
        prefix=lambda w, output: os.path.splitext(output.qopt)[0],
        extra=config["params"]["ngsadmix"]["extra"],
        reps=config["params"]["ngsadmix"]["reps"],
        minreps=config["params"]["ngsadmix"]["minreps"],
        thresh=config["params"]["ngsadmix"]["thresh"],
        conv=config["params"]["ngsadmix"]["conv"],
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt * 2880,
    script:
        "../scripts/ngsadmix.sh"


rule plot_admix:
    """
    Plot admixture proportions for a given K.
    """
    input:
        qopts=expand(
            "results/datasets/{{dataset}}/analyses/ngsadmix/{{dataset}}.{{ref}}_{{population}}{{dp}}_{{sites}}-filts_K{kvalue}.qopt",
            kvalue=config["params"]["ngsadmix"]["kvalues"],
        ),
        poplist="results/datasets/{dataset}/poplists/{dataset}_{population}{dp}.indiv.list",
        optwrap=expand(
            "results/datasets/{{dataset}}/analyses/ngsadmix/{{dataset}}.{{ref}}_{{population}}{{dp}}_{{sites}}-filts_K{kvalue}_optimization_wrapper.log",
            kvalue=config["params"]["ngsadmix"]["kvalues"],
        ),
    output:
        report(
            "results/datasets/{dataset}/plots/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts.pdf",
            category="03.2 Admixture",
            labels=lambda w: {
                "Filter": "{sites}",
                **dp_report(w),
                "Type": "Admixture plot",
            },
        ),
        convsumm="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts.convergence_summary.tsv",
    log:
        "logs/{dataset}/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_plot.log",
    benchmark:
        "benchmarks/{dataset}/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_plot.log"
    params:
        kvals=config["params"]["ngsadmix"]["kvalues"],
        thresh=config["params"]["ngsadmix"]["thresh"],
        conv=config["params"]["ngsadmix"]["conv"],
    container:
        r_container
    script:
        "../scripts/plot_admix.R"


rule admix_convergence_table:
    """Produce table convergence summary to add to report"""
    input:
        "results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts.convergence_summary.tsv",
    output:
        report(
            "results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts.convergence_summary.html",
            category="03.2 Admixture",
            labels=lambda w: {
                "Filter": "{sites}",
                **dp_report(w),
                "Type": "Convergence Table",
            },
        ),
    log:
        "logs/{dataset}/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_tsv2html.log",
    benchmark:
        "benchmarks/{dataset}/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_tsv2html.log"
    container:
        r_container
    shadow:
        "minimal"
    script:
        "../scripts/tsv2html.R"


rule evalAdmix:
    """
    Perform evalAdmix analysis to get residuals on best fitting replicate for a given K.
    """
    input:
        beagle=expand(
            "results/datasets/{{dataset}}/beagles/pruned/{{dataset}}.{{ref}}_{{population}}{{dp}}_{{sites}}-filts.pruned_maxkbdist-{maxkb}_minr2-{r2}.beagle.gz",
            maxkb=config["params"]["ngsld"]["max_kb_dist_pruning_dataset"],
            r2=config["params"]["ngsld"]["pruning_min-weight_dataset"],
        ),
        qopt="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.qopt",
        fopt="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.fopt.gz",
    output:
        "results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.corres",
    log:
        "logs/{dataset}/evaladmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.log",
    benchmark:
        "benchmarks/{dataset}/evaladmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.log"
    container:
        evaladmix_container
    shell:
        """
        evalAdmix -beagle {input.beagle} -fname {input.fopt} \
            -qname {input.qopt} -o {output} -P {threads} &> {log}
        """


rule plot_evalAdmix:
    """
    Plot evalAdmix residuals.
    """
    input:
        corres="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.corres",
        qopt="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.qopt",
        pops="results/datasets/{dataset}/poplists/{dataset}_{population}{dp}.indiv.list",
    output:
        report(
            "results/datasets/{dataset}/plots/evaladmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}_evaladmix.png",
            category="03.2 Admixture",
            subcategory="evalAdmix",
            labels=lambda w: {
                "Filter": "{sites}",
                **dp_report(w),
                "K-value": "{kvalue}",
                "Type": "Admix Corr. Residuals Plot",
            },
        ),
    log:
        "logs/{dataset}/evaladmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}_plot.log",
    benchmark:
        "benchmarks/{dataset}/evaladmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}_plot.log"
    container:
        evaladmix_container
    script:
        "../scripts/plot_evaladmix.R"
