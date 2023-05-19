# Rules for generating beagle (SNP genotype likelihoods) and maf (minor allele
# frequency) files per population


rule angsd_doGlf2:
    """
    Generates beagle and minor allele frequency files for a given population and genome 
    chunk. Calls SNPs from the whole dataset, and uses these same sites across all 
    population beagle files, even if a population is fixed for a certain allele.
    """
    input:
        glf=get_glf,
        fai="results/ref/{ref}/{ref}.fa.fai",
        sites=get_snpset,
    output:
        beagle=temp(
            "results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.beagle.gz"
        ),
        maf=temp(
            "results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.mafs.gz"
        ),
        arg="results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.arg",
    log:
        "logs/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    container:
        angsd_container
    params:
        popopts=get_popopts,
        nind=lambda w: len(get_samples_from_pop(w.population)),
        out=lambda w, output: os.path.splitext(output.arg)[0],
    threads: lambda wildcards, attempt: attempt
    resources:
        time=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        angsd -doGlf 2 -glf10_text {input.glf} {params.popopts} -doMaf 1 \
            -nThreads {threads} -sites {input.sites[0]} -fai {input.fai} \
            -nInd {params.nind} -out {params.out} &> {log}
        """


rule merge_beagle:
    """
    Merge beagle files across genome chunks per population.
    """
    input:
        lambda w: expand(
            "results/datasets/{{dataset}}/beagles/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.beagle.gz",
            chunk=chunklist,
        ),
    output:
        beagle="results/datasets/{dataset}/beagles/{dataset}.{ref}_{population}{dp}_{sites}-filts.beagle.gz",
    log:
        "logs/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_{sites}-filts_merge-beagle.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_{sites}-filts_merge-beagle.log"
    conda:
        "../envs/shell.yaml"
    resources:
        time=lambda wildcards, attempt: attempt * 60,
    shell:
        r"""
        (set +o pipefail;
        zcat {input} | head -n 1 | gzip > {output}

        for f in {input}; do
            zcat $f | tail -n +2 | gzip | cat >> {output.beagle}
        done) 2> {log}
        """


rule merge_maf:
    """
    Merges maf files across genome chunks per population.
    """
    input:
        lambda w: expand(
            "results/datasets/{{dataset}}/beagles/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.mafs.gz",
            chunk=chunklist,
        ),
    output:
        maf="results/datasets/{dataset}/mafs/{dataset}.{ref}_{population}{dp}_{sites}-filts.mafs.gz",
    log:
        "logs/{dataset}/angsd/doGlf2/{dataset}{ref}_{population}{dp}_{sites}-filts_merge-mafs.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doGlf2/{dataset}{ref}_{population}{dp}_{sites}-filts_merge-mafs.log"
    conda:
        "../envs/shell.yaml"
    resources:
        time=lambda wildcards, attempt: attempt * 60,
    shell:
        r"""
        (set +o pipefail;
        zcat {input} | head -n 1 | gzip > {output.maf}

        for f in {input}; do
            zcat $f | tail -n +2 | gzip | cat >> {output.maf}
        done) 2> {log}
        """


rule snpset:
    """
    Extracts SNP coordinates from maf file to get a list of variable sites in the 
    dataset.
    """
    input:
        "results/datasets/{dataset}/mafs/{dataset}.{ref}_all{dp}_{sites}-filts.mafs.gz",
    output:
        "results/datasets/{dataset}/filters/snps/{dataset}.{ref}{dp}_{sites}-filts_snps.sites",
    log:
        "logs/{dataset}/filters/snps/{dataset}.{ref}{dp}_{sites}-filts_snps.log",
    benchmark:
        "benchmarks/{dataset}/filters/snps/{dataset}.{ref}{dp}_{sites}-filts_snps.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        zcat {input} | tail -n +2 | cut -f1-4 > {output} 2> {log}
        """
