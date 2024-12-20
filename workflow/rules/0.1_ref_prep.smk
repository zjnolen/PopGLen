# Rules for preparing a reference genome for use with various tools. This includes
# linking it to an appropriate directory, indexing it, and dividing contigs up amongst
# 'chunks', which improves parallelization of analyses.


localrules:
    link_ref,
    link_anc_ref,
    samtools_faidx,
    ref_chunking,


rule link_ref:
    """Link reference genome to results directory"""
    input:
        config["reference"]["fasta"],
    output:
        "results/ref/{ref}/{ref}.fa",
    log:
        "logs/ref/link_ref/{ref}.log",
    container:
        shell_container
    resources:
        runtime="5m",
    shell:
        """
        ln -sr {input} {output} 2> {log}
        """


if config["ancestral"]:

    rule link_anc_ref:
        """Link ancestral reference genome to results directory if it exists"""
        input:
            config["ancestral"],
        output:
            "results/ref/{ref}/{ref}.anc.fa",
        log:
            "logs/ref/link_ref/{ref}.anc.log",
        container:
            shell_container
        resources:
            runtime="5m",
        shell:
            """
            ln -sr {input} {output} 2> {log}
            """


rule bwa_index:
    """Index reference genome for bwa (mapping)"""
    input:
        "results/ref/{ref}/{ref}.fa",
    output:
        idx=multiext(
            "results/ref/{ref}/{ref}.fa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "logs/ref/bwa_index/{ref}.log",
    resources:
        runtime="120m",
    benchmark:
        "benchmarks/ref/bwa_index/{ref}.log"
    wrapper:
        "v4.0.0/bio/bwa/index"


rule samtools_faidx:
    """
    Index reference genome using samtools (fai index used by several tools)
    """
    input:
        "results/ref/{ref}/{prefix}.fa",
    output:
        "results/ref/{ref}/{prefix}.fa.fai",
    container:
        samtools_container
    log:
        "logs/ref/samtools_faidx/{ref}/{prefix}.log",
    benchmark:
        "benchmarks/ref/samtools_faidx/{ref}/{prefix}.log"
    resources:
        runtime="10m",
    shell:
        """
        samtools faidx {input} 2> {log}
        """


rule ref_chunking:
    """Create regions files for running ANGSD in chunks along the genome"""
    input:
        "results/ref/{ref}/{ref}.fa",
    output:
        "results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
    log:
        "logs/{dataset}/ref/chunking/{ref}_chunk{chunk}.rf",
    container:
        shell_container
    params:
        contigs=lambda w: chunks[int(w.chunk) - 1].index.tolist(),
    resources:
        runtime="5m",
    shell:
        r"""
        echo {params.contigs} | tr " " "\n" > {output} 2> {log}
        """


rule picard_dict:
    """Create a dictionary for the reference using Picard for use with GATK"""
    input:
        "results/ref/{ref}/{ref}.fa",
    output:
        "results/ref/{ref}/{ref}.dict",
    log:
        "logs/ref/picard_dict/{ref}.log",
    benchmark:
        "benchmarks/ref/picard_dict/{ref}.log"
    resources:
        runtime="10m",
    wrapper:
        "v4.0.0/bio/picard/createsequencedictionary"
