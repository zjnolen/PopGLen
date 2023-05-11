# Estimates of inbreeding in the form of F_ROH - an inbreeding coefficient representing
# the proportion of the genome in runs of homozygosity greater than a certain length


rule ngsf_hmm:
    """
    Estimate IBD tracts within individual genomes.
    """
    input:
        beagle="results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_{population}{dp}_pruned.beagle.gz",
    output:
        ibd="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}.ibd",
        indF="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}.indF",
        pos="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}.pos",
    log:
        "logs/{dataset}/ngsF-HMM/{dataset}.{ref}_{population}{dp}.log",
    benchmark:
        "benchmarks/{dataset}/ngsF-HMM/{dataset}.{ref}_{population}{dp}.log"
    container:
        ngsf_hmm_container
    params:
        out=lambda w, output: os.path.splitext(output.pos)[0],
        nind=lambda w: len(get_samples_from_pop(w.population)),
    threads: lambda wildcards, attempt: attempt * 10
    resources:
        time=lambda wildcards, attempt: attempt * 2880,
    shell:
        r"""
        (zcat {input.beagle} | awk '{{print $1}}' | sed 's/\(.*\)_/\1\t/' \
            | tail -n +2 > {output.pos}
        
        nsites=$(cat {output.pos} | wc -l)

        nind={params.nind}

        export TMPDIR={resources.tmpdir}

        workflow/scripts/ngsF-HMM.sh --geno {input.beagle} --n_ind $nind \
            --n_sites $nsites --pos {output.pos} --lkl \
            --out {params.out} --n_threads {threads}) 2> {log}
        """


rule convert_ibd:
    """
    Converts ngsF-HMM ibd format into a list of runs of homozygosity.
    """
    input:
        ibd="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}.ibd",
        inds="results/datasets/{dataset}/poplists/{dataset}_{population}.indiv.list",
        pos="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}.pos",
    output:
        roh="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}.roh",
    log:
        "logs/{dataset}/ngsF-HMM/{dataset}.{ref}_{population}{dp}_convert_ibd.log",
    benchmark:
        "benchmarks/{dataset}/ngsF-HMM/{dataset}.{ref}_{population}{dp}_convert_ibd.log"
    container:
        ngsf_hmm_container
    shadow:
        "minimal"
    shell:
        """
        (tail -n +2 {input.inds} > {input.inds}.tmp

        # convert_ibd.pl drops warnings when chromosomes have non-numeric 
        # names, and may not handle them well. As a work around, a numeric 
        # index of chromosome names is created, then the names replaced with 
        # these numbers. Then returns their names after the conversion is 
        # complete.

        # first, create the index file
        n_contigs=$(awk '{{print $1}}' {input.pos} | uniq | wc -l) 2>> {log}
        awk '{{print $1}}' {input.pos} | uniq > {input.pos}.contigs 2>> {log}
        seq $n_contigs > {input.pos}.contigs.idx 2>> {log}
        
        # create index based pos file, adapted from:
        # Glenn Jackman - https://stackoverflow.com/a/7198895
        awk '
            FILENAME == ARGV[1] {{ listA[$1] = FNR; next }}
            FILENAME == ARGV[2] {{ listB[FNR] = $1; next }}
            {{
                for (i = 1; i <= NF; i++) {{
                    if ($1 in listA) {{
                           $1 = listB[listA[$i]]
                    }}
                }}
                print
            }}' \
        {input.pos}.contigs {input.pos}.contigs.idx {input.pos} \
            > {input.pos}.tmp

        convert_ibd.pl --pos_file {input.pos}.tmp \
            --ind_file {input.inds}.tmp    --ibd_pos_file {input.ibd} \
            > {output.roh}.tmp
        
        # Revert back to original contig name:
        awk '
            FILENAME == ARGV[1] {{ listA[$1] = FNR; next }}
            FILENAME == ARGV[2] {{ listB[FNR] = $1; next }}
            {{
                for (i = 1; i <= NF; i++) {{
                    if ($1 in listA) {{
                           $1 = listB[listA[$i]]
                    }}
                }}
                print
            }}' \
        {input.pos}.contigs.idx {input.pos}.contigs {output.roh}.tmp \
            > {output.roh}) 2> {log}
        """


rule plot_froh:
    """
    Plots average F_ROH per population for three IBD length thresholds (currently 
    hardcoded, will improve to accept options).
    """
    input:
        roh=expand(
            "results/datasets/{{dataset}}/analyses/ngsF-HMM/{{dataset}}.{{ref}}_{population}{{dp}}.roh",
            population=pop_list,
        ),
        inds="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
        autos=get_auto_sum,
    output:
        report(
            expand(
                "results/datasets/{{dataset}}/plots/inbreeding/{{dataset}}.{{ref}}_all{{dp}}.{stat}.pdf",
                stat=["froh", "rohreg"],
            ),
            category="Inbreeding",
        ),
    log:
        "logs/{dataset}/ngsF-HMM/{dataset}.{ref}_all{dp}_plot.log",
    benchmark:
        "benchmarks/{dataset}/ngsF-HMM/{dataset}.{ref}_all{dp}_plot.log"
    conda:
        "../envs/r.yaml"
    params:
        popnames=pop_list,
        outpre=lambda w, output: output[0].removesuffix(".froh.pdf"),
    script:
        "../scripts/plot_Froh.R"
