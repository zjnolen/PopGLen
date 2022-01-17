rule bwa_mem:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output
    output:
        temp(results+"/intermediate/mapped/{sample}.mem.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=genome_file(),
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    resources:
        time="24:00:00"
    wrapper:
        "0.84.0/bio/bwa/mem"

rule mark_duplicates:
    input:
        results+"/intermediate/mapped/{sample}.mem.bam"
    output:
        bam=protected(results+"/bam_dedup/{sample}.mem.bam"),
        metrics=results+"/mapping/dedup/{sample}.mem.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
#        config["params"]["picard"]["MarkDuplicates"]
    resources:
        time="12:00:00"
    wrapper:
        "0.84.0/bio/picard/markduplicates"

rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        "logs/samtools/index/{prefix}.bam"
    wrapper:
        "0.84.0/bio/samtools/index"