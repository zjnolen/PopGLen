localrules: link_ref, ref_chunking

rule link_ref:
    input:
        config["reference"]["fasta"]
    output:
        "results/ref/{ref}/{ref}.fa"
    log:
        "logs/ref/link_ref/{ref}.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        ln -sr {input} {output} 2> {log}
        """

rule bwa_index:
    input:
        "results/ref/{ref}/{ref}.fa"
    output:
        multiext("results/ref/{ref}/{ref}.fa",".amb",".ann",".bwt",".pac",".sa")
    log:
        "logs/ref/bwa_index/{ref}.log"
    conda:
        "../envs/mapping.yaml"
    resources:
        time=120
    shell:
        """
        bwa index {input} 2> {log}
        """

rule samtools_faidx:
    input:
        "results/ref/{ref}/{ref}.fa"
    output:
        "results/ref/{ref}/{ref}.fa.fai"
    log:
        "logs/ref/samtools_faidx/{ref}.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools faidx {input} 2> {log}
        """

rule ref_chunking:
    input:
        "results/ref/{ref}/{ref}.fa"
    output:
        "results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf"
    log:
        "logs/{dataset}/ref/chunking/{ref}_chunk{chunk}.rf"
    conda:
        "../envs/shell.yaml"
    params:
        contigs = lambda w: chunks[int(w.chunk) - 1].index.tolist()
    shell:
        r"""
        echo {params.contigs} | tr " " "\n" > {output} 2> {log}
        """

rule picard_dict:
    input:
        "results/ref/{ref}/{ref}.fa"
    output:
        "results/ref/{ref}/{ref}.dict"
    log:
        "logs/ref/picard_dict/{ref}.log"
    conda:
        "../envs/picard.yaml"
    shell:
        r"""
        picard CreateSequenceDictionary -Xmx{resources.mem_mb}m \
            R={input} O={output} 2> {log}
        """