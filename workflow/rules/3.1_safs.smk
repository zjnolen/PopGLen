# Rules for creating site allele frequency files used by ANGSD


rule angsd_doSaf:
    """
    Generate a site allele frequency file for a given population and genome chunk.
    """
    input:
        glf="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.glf.gz",
        fai="results/ref/{ref}/{ref}.fa.fai",
        anc="results/ref/{ref}/{ref}.fa",
    output:
        saf=temp(
            "results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.saf.gz"
        ),
        safidx=temp(
            "results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.saf.idx"
        ),
        safpos=temp(
            "results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.saf.pos.gz"
        ),
        arg="results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.arg",
    log:
        "logs/{dataset}/angsd/doSaf/{dataset}.{ref}_{population}{dp}_chunk{chunk}.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doSaf/{dataset}.{ref}_{population}{dp}_chunk{chunk}.log"
    container:
        angsd_container
    params:
        nind=lambda w: len(get_samples_from_pop(w.population)),
        out=lambda w, output: os.path.splitext(output.arg)[0],
    resources:
        time=lambda wildcards, attempt: attempt * 180,
    threads: lambda wildcards, attempt: attempt * 2
    shell:
        """
        angsd -doSaf 1 -glf10_text {input.glf} -anc {input.anc} \
            -nThreads 1 -fai {input.fai} -nInd {params.nind} \
            -out {params.out} &> {log}
        """


rule realSFS_catsaf:
    """
    Merge genome chunks into one file per population.
    """
    input:
        safs=lambda w: expand(
            "results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.saf.idx",
            chunk=chunklist,
        ),
        safgz=lambda w: expand(
            "results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.saf.gz",
            chunk=chunklist,
        ),
        safposgz=lambda w: expand(
            "results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.saf.pos.gz",
            chunk=chunklist,
        ),
    output:
        multiext(
            "results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}.saf",
            ".idx",
            ".pos.gz",
            ".gz",
        ),
    log:
        "logs/{dataset}/realSFS/cat/{dataset}.{ref}_{population}{dp}.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/cat/{dataset}.{ref}_{population}{dp}.log"
    container:
        angsd_container
    params:
        out=lambda w, output: output[0].removesuffix(".saf.idx"),
    resources:
        time=lambda wildcards, attempt: attempt * 60,
    shell:
        """
        realSFS cat {input.safs} -P 1 -outnames {params.out} 2> {log}
        """
