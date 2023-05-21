# Rules for creating site allele frequency files used by ANGSD


rule angsd_doSaf:
    """
    Generate a site allele frequency file for a given population and genome chunk.
    """
    input:
        glf=get_glf,
        fai="results/ref/{ref}/{ref}.fa.fai",
        anc="results/ref/{ref}/{ref}.fa",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.sites",
        idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.sites.idx",
    output:
        saf=temp(
            "results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.saf.gz"
        ),
        safidx=temp(
            "results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.saf.idx"
        ),
        safpos=temp(
            "results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.saf.pos.gz"
        ),
        arg="results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.arg",
    log:
        "logs/{dataset}/angsd/doSaf/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doSaf/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    container:
        angsd_container
    params:
        nind=lambda w: len(get_samples_from_pop(w.population)),
        out=lambda w, output: os.path.splitext(output.arg)[0],
    resources:
        runtime=lambda wildcards, attempt: attempt * 180,
    threads: lambda wildcards, attempt: attempt * 2
    shell:
        """
        angsd -doSaf 1 -glf10_text {input.glf} -anc {input.anc} -nThreads 1 \
            -fai {input.fai} -nInd {params.nind} -sites {input.sites} \
            -out {params.out} &> {log}
        """


rule realSFS_catsaf:
    """
    Merge genome chunks into one file per population.
    """
    input:
        safs=lambda w: expand(
            "results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.saf.idx",
            chunk=chunklist,
        ),
        safgz=lambda w: expand(
            "results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.saf.gz",
            chunk=chunklist,
        ),
        safposgz=lambda w: expand(
            "results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.saf.pos.gz",
            chunk=chunklist,
        ),
    output:
        multiext(
            "results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf",
            ".idx",
            ".pos.gz",
            ".gz",
        ),
    log:
        "logs/{dataset}/realSFS/cat/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/cat/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    container:
        angsd_container
    params:
        out=lambda w, output: output[0].removesuffix(".saf.idx"),
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        """
        realSFS cat {input.safs} -P 1 -outnames {params.out} 2> {log}
        """
