rule combine_LD_files:
    input:
        ldgz=expand(
            "results/datasets/{{dataset}}/analyses/ngsLD/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.ld.gz",
            chunk=chunklist,
        ),
    output:
        ldgz="results/datasets/{dataset}/analyses/ngsLD/{dataset}.{ref}_{population}{dp}_{sites}-filts.ld.gz",
    log:
        "logs/{dataset}/ngsLD/combine_LD_files/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/ngsLD/combine_LD_files/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        cat {input.ldgz} > {output.ldgz} 2> {log}
        """


rule fit_LD_decay:
    input:
        "results/datasets/{dataset}/analyses/ngsLD/{dataset}.{ref}_{population}{dp}_{sites}-filts.ld.gz"
    output:
        plot="results/datasets/{dataset}/plots/LD_decay/{dataset}.{ref}_{population}{dp}_{sites}-filts.LDdecay.svg",
        fit="results/datasets/{dataset}/analyses/ngsLD/{dataset}.{ref}_{population}{dp}_{sites}-filts.ld-decay-fit.txt"
    log:
        "logs/{dataset}/ngsLD/fit_LD_decay/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/ngsLD/fit_LD_decay/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    container:
        ngsld_container
    threads: lambda w, attempt: attempt
    resources:
        runtime=lambda w, attempt: attempt * 120,
    params:
        extra=config["params"]["ngsld"]["fit_extra"],
        nind=get_ngsld_n
    shell:
        """
        echo {input} | fit_LDdecay.R {params.extra} {params.nind} --out {output.plot} \
                > {output.fit} 2> {log}
        """
