# Rules for calculating genotype likelihoods from sequencing data using ANGSD


localrules:
    angsd_makeBamlist,
    popfile,


rule angsd_makeBamlist:
    """
    Create list of bam files to be input into calculation.
    """
    input:
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
    output:
        "results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
    log:
        "logs/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.log",
    conda:
        "../envs/shell.yaml"
    shell:
        """
        (readlink -f {input.bams} | perl -pe 'chomp if eof') > {output} 2> {log}
        """


rule popfile:
    """
    Create list of samples per population to preserve order of samples as sample names 
    are not propagated into downstream files.
    """
    output:
        inds="results/datasets/{dataset}/poplists/{dataset}_{population}.indiv.list",
    log:
        "logs/{dataset}/poplists/{dataset}_{population}_makelist.log",
    conda:
        "../envs/python.yaml"
    params:
        samplelist=samples,
        inds=lambda w: get_samples_from_pop(w.population),
    script:
        "../scripts/make_popfile.py"


rule angsd_sites_index:
    """
    Index sites files used for limiting calculations to specific positions.
    """
    input:
        "results/datasets/{dataset}/filters/{prefix}.sites",
    output:
        multiext("results/datasets/{dataset}/filters/{prefix}.sites", ".idx", ".bin"),
    log:
        "logs/{dataset}/angsd/sites_index/{prefix}.log",
    benchmark:
        "benchmarks/{dataset}/angsd/sites_index/{prefix}.log"
    container:
        angsd_container
    shell:
        """
        angsd sites index {input} 2> {log}
        """


rule angsd_doGlf4:
    """
    Calculate genotype likelihoods for filtered sites across all samples. Output to 
    text format for easy manipulation. Performed on subsets of the genome.
    """
    input:
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.sites",
        idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.sites.idx",
    output:
        glf="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.glf.gz",
        arg="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.arg",
    log:
        "logs/{dataset}/angsd/doGlf4/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doGlf4/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    wildcard_constraints:
        population="all",
    container:
        angsd_container
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: lambda wildcards, attempt: attempt * 2
    shell:
        """
        angsd -doGlf 4 -bam {input.bamlist} -GL {params.gl_model} \
            -ref {input.ref} -nThreads {threads} {params.extra} \
            -minMapQ {params.mapQ} -minQ {params.baseQ} -sites {input.sites} \
            -rf {input.regions} -out {params.out} &> {log}
        """


rule sampleglf:
    """
    Extract glf of single sample from dataset glf.
    """
    input:
        glf="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_all{dp}_chunk{chunk}_{sites}-filts.glf.gz",
    output:
        glf="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{sample}{dp}_chunk{chunk}_{sites}-filts.glf.gz",
    log:
        "logs/{dataset}/angsd/sampleglf/{dataset}.{ref}_{sample}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/sampleglf/{dataset}.{ref}_{sample}{dp}_chunk{chunk}_{sites}-filts.log"
    conda:
        "../envs/shell.yaml"
    params:
        start=lambda w: samples.index.values.tolist().index(w.sample),
        end=lambda w: str(samples.index.values.tolist().index(w.sample) + 1),
    shell:
        """
        zcat {input.glf} | cut -f1-2,{params.start}3-{params.end}2 | \
            gzip > {output} 2> {log}
        """


rule popglf:
    """
    Combine sample glfs into a population level glf.
    """
    input:
        sample_glfs=lambda w: expand(
            "results/datasets/{{dataset}}/glfs/chunks/{{dataset}}.{{ref}}_{sample}{{dp}}_chunk{{chunk}}_{{sites}}-filts.glf.gz",
            sample=get_samples_from_pop(w.population),
        ),
    output:
        glf="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.glf.gz",
    log:
        "logs/{dataset}/angsd/popglf/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/popglf/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    wildcard_constraints:
        population="|".join(
            [i for i in samples.index.tolist()]
            + [i for i in samples.population.values.tolist()]
            + [i for i in samples.depth.values.tolist()]
        ),
    conda:
        "../envs/shell.yaml"
    params:
        tmpfile="{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.glf",
    resources:
        runtime=360,
    shell:
        """
        (zcat {input.sample_glfs[0]} | cut -f1-2 > {resources.tmpdir}/{params.tmpfile}
        for i in {input.sample_glfs}; do
            echo "Adding $i to glf..."
            zcat $i | cut -f3-12 | paste -d '\t' {resources.tmpdir}/{params.tmpfile} - > {resources.tmpdir}/{params.tmpfile}.tmp
            mv {resources.tmpdir}/{params.tmpfile}.tmp {resources.tmpdir}/{params.tmpfile}
        done
        echo "Gzipping final glf..."
        gzip -c {resources.tmpdir}/{params.tmpfile} > {output.glf}) &> {log}
        """
