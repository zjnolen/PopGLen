# Estimates linkage disequilibrium between pairs of SNPs. If requested, saves a
# subsampling of these estimates (i.e. for estimating Ne using LD, investigating LD
# accross genome, etc.)


rule ngsLD_estLD:
    """
    Estimates pairwise linkage disequilibrium between SNPs.
    """
    input:
        beagle="results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.beagle.gz",
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
    output:
        ld=temp(
            "results/datasets/{dataset}/{path}/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.ld_maxkbdist-{maxkb}_rndsample-{rndsmp}.gz"
        ),
        pos=temp(
            "results/datasets/{dataset}/{path}/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.ld_maxkbdist-{maxkb}_rndsample-{rndsmp}.pos",
        ),
    log:
        "logs/{dataset}/ngsLD/estLD/{path}/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.ld_maxkbdist-{maxkb}_rndsample-{rndsmp}.log",
    benchmark:
        "benchmarks/{dataset}/ngsLD/estLD/{path}/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.ld_maxkbdist-{maxkb}_rndsample-{rndsmp}.log"
    container:
        ngsld_container
    threads: lambda wildcards, attempt: attempt
    wildcard_constraints:
        path="beagles/pruned/ngsLD|analyses/ngsLD/chunks",
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        r"""
        (zcat {input.beagle} | awk '{{print $1}}' | sed 's/\(.*\)_/\1\t/' \
            | tail -n +2 > {output.pos}
        
        nsites=$(cat {output.pos} | wc -l)

        nind=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
        if [[ $nsites == 0 ]]; then
            > {output.ld}
        else
            ngsLD --geno {input.beagle} --n_ind $nind --n_sites $nsites \
                --pos {output.pos} --probs --n_threads {threads} \
                --max_kb_dist {wildcards.maxkb} --rnd_sample {wildcards.rndsmp} \
                | gzip > {output.ld}
        fi) 2> {log}
        """


rule combine_LD_files:
    input:
        ldgz=expand(
            "results/datasets/{{dataset}}/analyses/ngsLD/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.ld_maxkbdist-{{maxkb}}_rndsample-{{rndsmp}}.gz",
            chunk=chunklist,
        ),
    output:
        ldgz="results/datasets/{dataset}/analyses/ngsLD/{dataset}.{ref}_{population}{dp}_{sites}-filts.ld_maxkbdist-{maxkb}_rndsample-{rndsmp}.gz",
    log:
        "logs/{dataset}/ngsLD/combine_LD_files/{dataset}.{ref}_{population}{dp}_{sites}-filts.ld_maxkbdist-{maxkb}_rndsample-{rndsmp}.log",
    benchmark:
        "benchmarks/{dataset}/ngsLD/combine_LD_files/{dataset}.{ref}_{population}{dp}_{sites}-filts.ld_maxkbdist-{maxkb}_rndsample-{rndsmp}.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        cat {input.ldgz} > {output.ldgz} 2> {log}
        """
