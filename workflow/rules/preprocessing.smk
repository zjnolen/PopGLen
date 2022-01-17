# adapter removal, quality trimming, and optionally pair collapsing with fastp

rule fastp_pe:
    input:
        sample=get_raw_fastq
    output:
        trimmed=expand("data/fastq_clean/{{sample}}.{read}.fastq.gz", read=["R1","R2"]),
        merged="data/fastq_clean/{sample}.merged.fastq.gz",
        html=results+"/preprocessing/fastp_reports/{sample}_paired.html",
        json=results+"/preprocessing/fastp_reports/{sample}_paired.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        extra="--merge"
    threads: 2
    resources:
        time="04:00:00"
    wrapper:
        "0.84.0/bio/fastp"