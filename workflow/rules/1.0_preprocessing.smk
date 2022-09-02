# adapter removal, quality trimming, and optionally pair collapsing with fastp

def get_raw_fastq(wildcards):
    # Checks if raw sequencing data is specified locally and sets fastp 
    # input as either raw local files or to download from SRA
    unit = units.loc[wildcards.sample, ["fq1", "fq2"]]
    return [unit.fq1,unit.fq2]

rule fastp_mergedout:
    input:
        sample=get_raw_fastq
    output:
        merged=results+"/preprocessing/fastp/{sample}.merged.fastq.gz",
        html=results+"/preprocessing/fastp/{sample}_paired.html",
        json=results+"/preprocessing/fastp/{sample}_paired.json"
    log:
        logs+"/fastp/{sample}.merged.log"
    conda:
        "../envs/fastp.yaml"
    params:
        extra="-p -g --overlap_len_require 15"
    threads: lambda wildcards, attempt: attempt*2
    resources:
        time=lambda wildcards, attempt: attempt*240
    shell:
        """
        fastp -w {threads} {params.extra} -m -i {input[0]} \
            -I {input[1]} --merged_out {output.merged} \
            -j {output.json} -h {output.html}
        """

rule fastp_pairedout:
    input:
        sample=get_raw_fastq
    output:
        paired=expand(results+"/preprocessing/fastp/{{sample}}.{read}.fastq.gz",
            read=["R1","R2"]),
        html=results+"/preprocessing/fastp/{sample}_paired.html",
        json=results+"/preprocessing/fastp/{sample}_paired.json"
    log:
        logs+"/fastp/{sample}.paired.log"
    conda:
        "../envs/fastp.yaml"
    params:
        extra="-p -g --overlap_len_require 15"
    threads: lambda wildcards, attempt: attempt*2
    resources:
        time=lambda wildcards, attempt: attempt*240
    shell:
        """
        fastp -w {threads} {params.extra} -i {input[0]} \
            -I {input[1]} -o {output.paired[0]} -O {output.paired[1]} \
            -j {output.json} -h {output.html}
        """

# rule fastp_pe:
#     input:
#         sample=get_raw_fastq
#     output:
#         trimmed=expand(results+"/preprocessing/fastp/{{sample}}.{read}.fastq.gz",
#             read=["R1","R2"]),
#         merged=results+"/preprocessing/fastp/{sample}.merged.fastq.gz",
#         unpaired=results+"/preprocessing/fastp/{sample}.singletons.fastq.gz",
#         html=results+"/preprocessing/fastp/{sample}_paired.html",
#         json=results+"/preprocessing/fastp/{sample}_paired.json"
#     log:
#         logs+"/fastp/{sample}.log"
#     params:
#         extra="--merge -p -g"
#     threads: lambda wildcards, attempt: attempt*2
#     resources:
#         time=lambda wildcards, attempt: attempt*240
#     wrapper:
#         "0.84.0/bio/fastp"