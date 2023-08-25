# Rules for preparing a reference genome for use with various tools. This includes
# linking it to an appropriate directory, indexing it, and dividing contigs up amongst
# 'chunks', which improves parallelization of analyses.


localrules:
    link_ref,
    ref_chunking,


rule link_ref:
    """Link reference genome to results directory"""
    input:
        config["reference"]["fasta"],
    output:
        "results/ref/{ref}/{ref}.fa",
    log:
        "logs/ref/link_ref/{ref}.log",
    conda:
        "../envs/shell.yaml"
    shell:
        """
        ln -sr {input} {output} 2> {log}
        """


rule bwa_index:
    """Index reference genome for bwa (mapping)"""
    input:
        "results/ref/{ref}/{ref}.fa",
    output:
        multiext("results/ref/{ref}/{ref}.fa", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    log:
        "logs/ref/bwa_index/{ref}.log",
    resources:
        runtime=120,
    benchmark:
        "benchmarks/ref/bwa_index/{ref}.log"
    wrapper:
        "v2.6.0/bio/bwa-mem2/index"


rule samtools_faidx:
    """Index reference genome using samtools (fai index used by several tools)"""
    input:
        "results/ref/{ref}/{ref}.fa",
    output:
        "results/ref/{ref}/{ref}.fa.fai",
    log:
        "logs/ref/samtools_faidx/{ref}.log",
    benchmark:
        "benchmarks/ref/samtools_faidx/{ref}.log"
    wrapper:
        "v2.4.0/bio/samtools/faidx"


rule ref_chunking:
    """Create regions files for running ANGSD in chunks along the genome"""
    input:
        "results/ref/{ref}/{ref}.fa",
    output:
        "results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
    log:
        "logs/{dataset}/ref/chunking/{ref}_chunk{chunk}.rf",
    conda:
        "../envs/shell.yaml"
    params:
        contigs=lambda w: chunks[int(w.chunk) - 1].index.tolist(),
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
    wrapper:
        "v2.4.0/bio/picard/createsequencedictionary"
