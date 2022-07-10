rule bwa_mem:
    input:
        merged=rules.fastp_pe.output.merged,
        unpaired=rules.fastp_pe.output.unpaired,
        paired=rules.fastp_pe.output.trimmed,
        idx=rules.bwa_index.output
    output:
        SEfastq=results+"/mapping/mapped/{sample}.SE.fastq.gz",
        singlebam=results+"/mapping/mapped/{sample}.singles.bam",
        pairbam=results+"/mapping/mapped/{sample}.pairs.bam",
        allbam=results+"/mapping/mapped/{sample}.bam"
    log:
        logs + "/bwa_mem/{sample}.log"
    params:
        index=REF,
        rg=get_read_group
    conda:
        "../envs/mapping.yaml"
    threads: lambda wildcards, attempt: attempt*4
    resources:
        time=lambda wildcards, attempt: attempt*720
    shell:
        """
        # Maps merged + unpaired (SE) and trimmed (PE) reads separately, then 
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
            {output.SEfastq} | \
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
        results+"/mapping/mapped/{sample}.bam"
    output:
        bam=protected(results+"/mapping/dedup/{sample}.bam"),
        metrics=results+"/mapping/dedup/{sample}.metrics.txt"
    log:
        logs + "/picard/dedup/{sample}.log"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"]
    threads: lambda wildcards, attempt: attempt*2
    resources:
        time=lambda wildcards, attempt: attempt*720
    wrapper:
        "0.84.0/bio/picard/markduplicates"

rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        "0.84.0/bio/samtools/index"

rule samtools_subsample:
    input:
        bam=results+"/mapping/dedup/{sample}.bam"
    output:
        results+"/mapping/dedup/{sample}{dp}.bam"
    log:
        logs+"/samtools/subsample/{sample}{dp}.log"
    conda:
        "../envs/samtools.yaml"
    resources:
        time="06:00:00"
    shell:
        """
        dp=$(samtools depth -a {input.bam} | awk '{{sum+=$3}} END \
            {{print sum/NR}}')
        
        subdp=$(echo {wildcards.dp} | sed 's/.dp//g')

        prop=$(echo "$dp $subdp" | awk '{{print 1 / ($1 / $2)}}')

        if [ `awk 'BEGIN {{print ('$prop' <= 1.0)}}'` = 1 ]; then
            propdec=$(echo $prop | awk -F "." '{{print $2}}')
            samtools view -h -s ${{RANDOM}}.${{propdec}} -@ {threads} \
                -b {input.bam} > {output} 2> {log}
        else
            original=$(readlink -f {input.bam})
            ln -sf $original {output}
        fi

        echo "Subsampled average depth:" >> {log}

        samtools depth -a {output} | awk '{{sum+=$3}} END \
            {{print sum/NR}}' &>> {log}
        """

# rule samtools_subsample:
#     input:
#         bam=results + "/mapping/dedup/{sample}.bam",
#         DP=results + "/depth/{sample}.depthMean"
#     output:
#         results + "/dedup/{sample}.dp{cov}.bam"
#     log:
#         logs + "samtools/subsample/{sample}.dp{cov}.log"
#     conda:
#         "../envs/samtools.yaml"
#     resources:
#         time="06:00:00"
#     shell:
#         """
#         cov=$(awk '{{print $1}}' {input.DP})
        
#         propcov=$(echo "$cov {wildcards.cov}" | awk '{{print 1 / ($1 / $2)}}')
        
#         if (( $(echo "$propcov > 1" | bc -l) )); then
#             ln {wildcards.sample}.bam {output}
#         else
#             samtools view -s $propcov -b {input.bam} > {output}
#         fi
#         """