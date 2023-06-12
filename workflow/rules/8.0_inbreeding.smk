# Estimates of inbreeding in the form of F_ROH - an inbreeding coefficient representing
# the proportion of the genome in runs of homozygosity greater than a certain length


rule ngsf_hmm:
    """
    Estimate IBD tracts within individual genomes.
    """
    input:
        beagle="results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_{population}{dp}_{sites}-filts_pruned.beagle.gz",
    output:
        ibd="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts.ibd",
        indF="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts.indF",
        pos="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts.pos",
    log:
        "logs/{dataset}/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    container:
        ngsf_hmm_container
    params:
        out=lambda w, output: os.path.splitext(output.pos)[0],
        nind=lambda w: len(get_samples_from_pop(w.population)),
    threads: lambda wildcards, attempt: attempt * 10
    resources:
        runtime=lambda wildcards, attempt: attempt * 2880,
    script:
        "../scripts/ngsF-HMM.sh"


rule convert_ibd:
    """
    Converts ngsF-HMM ibd format into a list of runs of homozygosity.
    """
    input:
        ibd="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts.ibd",
        inds="results/datasets/{dataset}/poplists/{dataset}_{population}.indiv.list",
        pos="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts.pos",
    output:
        roh="results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts.roh",
    log:
        "logs/{dataset}/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts_convert_ibd.log",
    benchmark:
        "benchmarks/{dataset}/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts_convert_ibd.log"
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
        n_contigs=$(awk '{{print $1}}' {input.pos} | uniq | wc -l)
        awk '{{print $1}}' {input.pos} | uniq > {input.pos}.contigs
        seq $n_contigs > {input.pos}.contigs.idx
        paste -d "\t" {input.pos}.contigs {input.pos}.contigs.idx > {input.pos}.index
        
        # create index based pos file, adapted from:
        # Jaypal Singh - https://stackoverflow.com/a/22253586
        awk 'NR==FNR{{a[$1]=$2;next}}{{$1=a[$1]}}1' {input.pos}.index {input.pos} \
            > {input.pos}.tmp

        convert_ibd.pl --pos_file {input.pos}.tmp \
            --ind_file {input.inds}.tmp    --ibd_pos_file {input.ibd} \
            > {output.roh}.tmp
        
        # Revert back to original contig name:
        awk 'NR==FNR{{a[$2]=$1;next}}{{$1=a[$1]}}1' {input.pos}.index {output.roh}.tmp \
            > {output.roh}) 2> {log}
        """


rule plot_froh:
    """
    Plots average F_ROH per population for three IBD length thresholds (currently 
    hardcoded, will improve to accept options).
    """
    input:
        roh=expand(
            "results/datasets/{{dataset}}/analyses/ngsF-HMM/{{dataset}}.{{ref}}_{population}{{dp}}_{{sites}}-filts.roh",
            population=pop_list,
        ),
        inds="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
        autos=get_auto_sum,
    output:
        plot=report(
            "results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.froh.pdf",
            category="Inbreeding",
            labels=lambda w: {"Filter": "{sites}", **dp_report(w), "Type": "Barplot"},
        ),
        tsv="results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.froh.tsv",
    log:
        "logs/{dataset}/ngsF-HMM/{dataset}.{ref}_all{dp}_{sites}-filts_plot.log",
    benchmark:
        "benchmarks/{dataset}/ngsF-HMM/{dataset}.{ref}_all{dp}_{sites}-filts_plot.log"
    conda:
        "../envs/r.yaml"
    params:
        popnames=pop_list,
        outpre=lambda w, output: output["plot"].removesuffix(".froh.pdf"),
    script:
        "../scripts/plot_Froh.R"
