# adapter removal, quality trimming, and optionally pair collapsing with fastp

rule fastp_pe:
    input:
        sample=get_raw_fastq
    output:
        trimmed=temp(expand(intermediate+"/fastp/{{sample}}.{read}.fastq.gz", read=["R1","R2"])),
        merged=temp(intermediate+"/fastp/{sample}.merged.fastq.gz"),
        html=results+"/fastp_reports/{sample}_paired.html",
        json=results+"/fastp_reports/{sample}_paired.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        extra="--merge"
    threads: 2
    resources:
        time="04:00:00"
    wrapper:
        "0.84.0/bio/fastp"