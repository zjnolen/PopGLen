# adapter removal, quality trimming, and optionally pair collapsing with fastp

rule fastp_pe:
    input:
        sample=get_raw_fastq
    output:
        trimmed=expand(intermediate+"/fastp/{{sample}}.{read}.fastq.gz",
            read=["R1","R2"]),
        merged=intermediate+"/fastp/{sample}.merged.fastq.gz",
        unpaired=intermediate+"/fastp/{sample}.singletons.fastq.gz",
        html=results+"/fastp/{sample}_paired.html",
        json=results+"/fastp/{sample}_paired.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        extra="--merge -p -g"
    threads: 2
    resources:
        time="04:00:00"
    wrapper:
        "0.84.0/bio/fastp"