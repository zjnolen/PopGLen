# Rules for adapter removal, quality trimming, and pair collapsing (for historical
# samples only) with fastp


rule ncbi_download:
    """
    If SRA accessions are the only input in a row of the units list, download
    the FASTQ files from NCBI.
    """
    output:
        temp("results/downloaded_fastq/{accession}_1.fastq.gz"),
        temp("results/downloaded_fastq/{accession}_2.fastq.gz"),
    log:
        "logs/download_fastq/{accession}.log",
    params:
        extra="--skip-technical -x",
    threads: 6
    resources:
        runtime="6h",
    wrapper:
        "v4.0.0/bio/sra-tools/fasterq-dump"


rule fastp_mergedout:
    """Process reads with fastp, collapsing overlapping read pairs"""
    input:
        unpack(get_raw_fastq),
    output:
        trimmed=temp(
            expand(
                "results/preprocessing/fastp/{{sample}}_{{unit}}_{{lib}}.{read}.uncollapsed.fastq.gz",
                read=["R1", "R2"],
            )
        ),
        merged=temp("results/preprocessing/fastp/{sample}_{unit}_{lib}.merged.fastq.gz"),
        html=report(
            "results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_merged.html",
            category="00 Quality Control",
            subcategory="1 Trimming Reports",
            labels={
                "Sample": "{sample}",
                "Unit": "{unit}",
                "Lib": "{lib}",
                "Type": "fastp Report",
            },
        ),
        json="results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_merged.json",
    log:
        "logs/preprocessing/fastp/{sample}_{unit}_{lib}.merged.log",
    benchmark:
        "benchmarks/preprocessing/fastp/{sample}_{unit}_{lib}.merged.log"
    params:
        extra=lambda w: config["params"]["fastp"]["extra"]
        + f" --merge --overlap_len_require {config['params']['fastp']['min_overlap_hist']}",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 480,
    wrapper:
        "v4.0.0/bio/fastp"


rule fastp_pairedout:
    """Process reads with fastp, don't collapse overlapping read pairs"""
    input:
        unpack(get_raw_fastq),
    output:
        trimmed=temp(
            expand(
                "results/preprocessing/fastp/{{sample}}_{{unit}}_{{lib}}.{read}.paired.fastq.gz",
                read=["R1", "R2"],
            )
        ),
        html=report(
            "results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_paired.html",
            category="00 Quality Control",
            subcategory="1 Trimming Reports",
            labels={
                "Sample": "{sample}",
                "Unit": "{unit}",
                "Lib": "{lib}",
                "Type": "fastp Report",
            },
        ),
        json="results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_paired.json",
    log:
        "logs/preprocessing/fastp/{sample}_{unit}_{lib}.paired.log",
    benchmark:
        "benchmarks/preprocessing/fastp/{sample}_{unit}_{lib}.paired.log"
    params:
        extra=lambda w: config["params"]["fastp"]["extra"],
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 480,
    wrapper:
        "v4.0.0/bio/fastp"


rule fastp_multiqc:
    """Combine fastp reports into a single MultiQC for all samples"""
    input:
        multiqc_input_fastp,
    output:
        report(
            "results/datasets/{dataset}/qc/fastp-trimming/fastp_all.{ref}_mqc.html",
            category="00 Quality Control",
            subcategory="1 Trimming Reports",
            labels={"Type": "MultiQC Report"},
        ),
    log:
        "logs/preprocessing/fastp/{dataset}.{ref}_mqc.log",
    container:
        multiqc_container
    resources:
        runtime="1h",
    shell:
        """
        multiqc --no-data-dir --filename {output} {input} 2> {log}
        """
