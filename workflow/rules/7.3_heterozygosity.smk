# Estimates of individual genome wide heterozygosity from site frequency spectrum


rule heterozygosity:
    """
    Calculates individual heterozygosity for all samples from single sample 1D SFS.
    Plots per population distributions of individual heterozygosity.
    """
    input:
        sfs=expand(
            "results/datasets/{{dataset}}/analyses/sfs/{{dataset}}.{{ref}}_{sample}{{dp}}_{{sites}}-filts.sfs",
            sample=samples.index,
        ),
        popfile="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
    output:
        "results/datasets/{dataset}/analyses/heterozygosity/{dataset}.{ref}_all{dp}_{sites}-filts_heterozygosity.tsv",
        report(
            "results/datasets/{dataset}/plots/heterozygosity/{dataset}.{ref}_all{dp}_{sites}-filts_heterozygosity.pdf",
            category="Heterozygosity",
            labels=lambda w: {"Filter": "{sites}", **dp_report(w), "Type": "Boxplot"},
        ),
    log:
        "logs/{dataset}/heterozygosity/{dataset}.{ref}_all{dp}_{sites}-filts_calc-plot.log",
    benchmark:
        "benchmarks/{dataset}/heterozygosity/{dataset}.{ref}_all{dp}_{sites}-filts_calc-plot.log"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/plot_heterozygosity.R"


rule heterozygosity_table:
    """
    Converts individual heterozygosity table to html for report.
    """
    input:
        "results/datasets/{dataset}/analyses/heterozygosity/{dataset}.{ref}_all{dp}_{sites}-filts_heterozygosity.tsv",
    output:
        report(
            "results/datasets/{dataset}/analyses/heterozygosity/{dataset}.{ref}_all{dp}_{sites}-filts_heterozygosity.html",
            category="Heterozygosity",
            labels=lambda w: {"Filter": "{sites}", **dp_report(w), "Type": "Table"},
        ),
    log:
        "logs/{dataset}/heterozygosity/{dataset}.{ref}_all{dp}_{sites}-filts_tsv2html.log",
    benchmark:
        "benchmarks/{dataset}/heterozygosity/{dataset}.{ref}_all{dp}_{sites}-filts_tsv2html.log"
    conda:
        "../envs/r-rectable.yaml"
    script:
        "../scripts/tsv2html.Rmd"
