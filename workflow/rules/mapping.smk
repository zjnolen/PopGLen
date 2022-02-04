rule bwa_mem:
    input:
        merged=rules.fastp_pe.output.merged,
        unpaired=rules.fastp_pe.output.unpaired,
        paired=rules.fastp_pe.output.trimmed,
        idx=rules.bwa_index.output
    output:
        SEfastq=intermediate+"/mapping/{sample}.SE.fastq.gz",
        singlebam=intermediate+"/mapping/{sample}.singles.mem.bam",
        pairbam=intermediate+"/mapping/{sample}.pairs.mem.bam",
        allbam=intermediate+"/mapping/{sample}.mem.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=genome_file(),
        rg=get_read_group
    conda:
        "../envs/mapping.yaml"
    threads: 4
    resources:
        time="24:00:00"
    shell:
        """
        # Maps merged, unpaired (SE) and trimmed (PE) reads separately, then 
        # merges them into a single bam file. Will maybe be updated to do all 
        # in parallel, and to be able to request only certain read types to be
        # used.

        # Combine merged and unpaired reads (all single ended now) to map them 
        # in one step.
        cat {input.merged} {input.unpaired} > {output.SEfastq} 2> {log}

        # Map SE reads
        bwa mem \
            -t {threads} \
            {params.rg} \
            {params.index} \
            {resources.tmpdir}/{wildcards.sample}.SE.fastq.gz | \
        samtools sort -o {output.singlebam} 2>> {log}

        # Map paired reads
        bwa mem \
            -t {threads} \
            {params.rg} \
            {params.index} \
            {input.paired} | \
        samtools sort -o {output.pairbam} 2>> {log}

        # Merge bam files for final output
        samtools merge \
            -@ {threads} \
            -o {output.allbam} \
            {output.singlebam} {output.pairbam} 2>> {log}
        """

rule mark_duplicates:
    input:
        intermediate+"/mapping/{sample}.mem.bam"
    output:
        bam=protected(results+"/dedup/{sample}.mem.bam"),
        metrics=results+"/dedup/{sample}.mem.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"]
    threads: 2
    resources:
        time="12:00:00",
        mem_mb=10240
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