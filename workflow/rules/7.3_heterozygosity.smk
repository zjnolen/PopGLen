# Estimates of individual genome wide heterozygosity from site frequency spectrum


rule heterozygosity:
    """
    Calculates individual heterozygosity for all samples from single sample 1D SFS.
    Plots per population distributions of individual heterozygosity.
    """
    input:
        sfs=expand(
            "results/datasets/{{dataset}}/analyses/sfs/{{dataset}}.{{ref}}_{sample}{{dp}}.sfs",
            sample=samples.index,
        ),
        popfile="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
    output:
        "results/datasets/{dataset}/analyses/heterozygosity/{dataset}.{ref}_all{dp}_heterozygosity.tsv",
        report(
            "results/datasets/{dataset}/plots/heterozygosity/{dataset}.{ref}_all{dp}_heterozygosity.pdf",
            category="Heterozygosity",
            labels={
                "Topic": "Heterozygosity",
                "Subsampling": "{dp}",
                "Type": "boxplot",
            },
        ),
    log:
        "logs/{dataset}/heterozygosity/{dataset}.{ref}_all{dp}_calc-plot.log",
    benchmark:
        "benchmarks/{dataset}/heterozygosity/{dataset}.{ref}_all{dp}_calc-plot.log"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/plot_heterozygosity.R"


rule heterozygosity_table:
    """
    Converts individual heterozygosity table to html for report.
    """
    input:
        "results/datasets/{dataset}/analyses/heterozygosity/{dataset}.{ref}_all{dp}_heterozygosity.tsv",
    output:
        report(
            "results/datasets/{dataset}/analyses/heterozygosity/{dataset}.{ref}_all{dp}_heterozygosity.html",
            category="Heterozygosity",
            labels={"Topic": "Heterozygosity", "Subsampling": "{dp}", "Type": "Table"},
        ),
    log:
        "logs/{dataset}/heterozygosity/{dataset}.{ref}_all{dp}_tsv2html.log",
    benchmark:
        "benchmarks/{dataset}/heterozygosity/{dataset}.{ref}_all{dp}_tsv2html.log"
    conda:
        "../envs/r-rectable.yaml"
    script:
        "../scripts/tsv2html.Rmd"
