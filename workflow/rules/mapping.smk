rule bwa_mem:
    input:
        reads=get_fastp_reads,
        idx=rules.bwa_index.output
    output:
        mergebam=temp(intermediate+"/mapped/{sample}.merge.mem.bam"),
        trimbam=temp(intermediate+"/mapped/{sample}.trim.mem.bam"),
        allbam=temp(intermediate+"/mapped/{sample}.mem.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=genome_file(),
        rg=get_read_group
    conda:
        "envs/mapping.yaml"
    threads: 4
    resources:
        time="24:00:00"
    shell:
        """
        # Maps merged (SE) and trimmed (PE) reads separately, then merges them
        # into a single bam file. Will later be updated to do all three stages
        # in parallel, and to be able to request only merged or trimmed reads
        $ to be used.

        # map merged reads
            
        bwa mem \
            -t {threads} \
            {params.rg} \
            {params.index} \
            {input.reads[0]} | \
        \
        samtools sort -o {output.mergebam}

        # map paired reads
        bwa mem \
            -t {threads} \
            {params.rg} \
            {params.index} \
            {input.reads[1]} {input.reads[2]} | \
        \
        samtools sort -o {output.trimbam}

        # merge bam files
        samtools merge \
            -@ {threads} \
            -o {output.allbam} \
            {output.mergebam} {output.trimbam}
        """

rule mark_duplicates:
    input:
        intermediate+"/mapped/{sample}.mem.bam"
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