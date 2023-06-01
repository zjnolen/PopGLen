# Rules for adapter removal, quality trimming, and pair collapsing (for historical
# samples only) with fastp


rule fastp_mergedout:
    """Process historical reads with fastp, collapsing overlapping read pairs"""
    input:
        unpack(get_raw_fastq),
    output:
        merged="results/preprocessing/fastp/{sample}.merged.fastq.gz",
        html=report(
            "results/preprocessing/qc/fastp/{sample}_paired.html",
            category="Quality Control",
            subcategory="Trimming Reports",
            labels={"Sample": "{sample}", "Type": "fastp Report"},
        ),
        json="results/preprocessing/qc/fastp/{sample}_paired.json",
    log:
        "logs/preprocessing/fastp/{sample}.merged.log",
    benchmark:
        "benchmarks/preprocessing/fastp/{sample}.merged.log"
    conda:
        "../envs/fastp.yaml"
    params:
        extra=config["params"]["fastp"]["extra"],
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 240,
    shell:
        """
        fastp -w {threads} {params.extra} -m -i {input.r1} \
            -I {input.r2} --merged_out {output.merged} \
            -j {output.json} -h {output.html} 2> {log}
        """


rule fastp_pairedout:
    """Process modern reads with fastp, trimming adapters and low quality bases"""
    input:
        unpack(get_raw_fastq),
    output:
        paired=expand(
            "results/preprocessing/fastp/{{sample}}.{read}.fastq.gz", read=["R1", "R2"]
        ),
        html=report(
            "results/preprocessing/qc/fastp/{sample}_paired.html",
            category="Quality Control",
            subcategory="Trimming Reports",
            labels={"Sample": "{sample}", "Type": "fastp Report"},
        ),
        json="results/preprocessing/qc/fastp/{sample}_paired.json",
    log:
        "logs/preprocessing/fastp/{sample}.paired.log",
    benchmark:
        "benchmarks/preprocessing/fastp/{sample}.paired.log"
    conda:
        "../envs/fastp.yaml"
    params:
        extra=config["params"]["fastp"]["extra"],
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 240,
    shell:
        """
        fastp -w {threads} {params.extra} -i {input.r1} \
            -I {input.r2} -o {output.paired[0]} -O {output.paired[1]} \
            -j {output.json} -h {output.html} 2> {log}
        """
