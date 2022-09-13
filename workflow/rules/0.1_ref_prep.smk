localrules: ref_chunking

rule bwa_index:
    input:
        REF
    output:
        multiext(REF,".amb",".ann",".bwt",".pac",".sa")
    conda:
        "../envs/mapping.yaml"
    params:
        extra=""
    resources:
        time=120
    shell:
        """
        bwa index {params.extra} {input}
        """

rule samtools_faidx:
    input:
        REF
    output:
        REF + ".fai"
    conda:
        "../envs/samtools.yaml"
    params:
        extra=""
    shell:
        """
        samtools faidx {params.extra} {input}
        """

rule ref_chunking:
    input:
        REF
    output:
        REF_DIR+"/beds/chunk{chunk}_"+str(config["chunk_size"])+"bp.rf"
    params:
        contigs = lambda w: chunks[int(w.chunk) - 1].index.tolist()
    shell:
        r"""echo {params.contigs} | tr " " "\n" > {output}"""

rule picard_dict:
    input:
        REF
    output:
        REF+".dict"
    conda:
        "../envs/picard.yaml"
    shell:
        r"""
        picard CreateSequenceDictionary -Xmx{resources.mem_mb}m \
            R={input} O={output}
        """