# rule bwa_mem:
#     input:
#         merged=rules.fastp_pe.output.merged,
#         unpaired=rules.fastp_pe.output.unpaired,
#         paired=rules.fastp_pe.output.trimmed,
#         idx=rules.bwa_index.output
#     output:
#         SEfastq=results+"/mapping/mapped/{sample}.SE.fastq.gz",
#         singlebam=results+"/mapping/mapped/{sample}.singles.bam",
#         pairbam=results+"/mapping/mapped/{sample}.pairs.bam",
#         allbam=results+"/mapping/mapped/{sample}.bam"
#     log:
#         logs + "/bwa_mem/{sample}.log"
#     params:
#         index=REF,
#         rg=get_read_group
#     conda:
#         "../envs/mapping.yaml"
#     threads: lambda wildcards, attempt: attempt*4
#     resources:
#         time=lambda wildcards, attempt: attempt*720
#     shell:
#         """
#         # Maps merged + unpaired (SE) and trimmed (PE) reads separately, then 
#         # merges them into a single bam file. Will maybe be updated to do all 
#         # in parallel, and to be able to request only certain read types to be
#         # used.

#         # Combine merged and unpaired reads (all single ended now) to map them 
#         # in one step.
#         cat {input.merged} {input.unpaired} > {output.SEfastq} 2> {log}

#         # Map SE reads
#         bwa mem \
#             -t {threads} \
#             {params.rg} \
#             {params.index} \
#             {output.SEfastq} | \
#         samtools sort -o {output.singlebam} 2>> {log}

#         # Map paired reads
#         bwa mem \
#             -t {threads} \
#             {params.rg} \
#             {params.index} \
#             {input.paired} | \
#         samtools sort -o {output.pairbam} 2>> {log}

#         # Merge bam files for final output
#         samtools merge \
#             -@ {threads} \
#             -o {output.allbam} \
#             {output.singlebam} {output.pairbam} 2>> {log}
#         """

rule bwa_mem_merged:
    input:
        merged=rules.fastp_mergedout.output.merged,
        idx=rules.bwa_index.output
    output:
        bam=results+"/mapping/mapped/{sample}.merged.bam"
    log:
        logs + "/bwa_mem/{sample}.merged.log"
    params:
        index=REF,
        rg=get_read_group
    conda:
        "../envs/mapping.yaml"
    shadow: "copy-minimal"
    threads: lambda wildcards, attempt: attempt*8
    resources:
        time=lambda wildcards, attempt: attempt*10080
    shell:
        """
        bwa mem \
            -t {threads} \
            {params.rg} \
            {params.index} \
            {input.merged} | \
        samtools sort -o {output.bam} 2> {log}
        """

rule bwa_mem_paired:
    input:
        paired=rules.fastp_pairedout.output.paired,
        idx=rules.bwa_index.output
    output:
        bam=results+"/mapping/mapped/{sample}.paired.bam"
    log:
        logs + "/bwa_mem/{sample}.paired.log"
    params:
        index=REF,
        rg=get_read_group
    conda:
        "../envs/mapping.yaml"
    shadow: "copy-minimal"
    threads: lambda wildcards, attempt: attempt*20
    resources:
        time=lambda wildcards, attempt: attempt*10080
    shell:
        """
        bwa mem \
            -t {threads} \
            {params.rg} \
            {params.index} \
            {input.paired} | \
        samtools sort -o {output.bam} 2> {log}
        """

rule mark_duplicates:
    input:
        results+"/mapping/mapped/{sample}.paired.bam"
    output:
        bam=results+"/mapping/dedup/{sample}.paired.rmdup.bam",
        metrics=results+"/mapping/dedup/{sample}.metrics.txt"
    log:
        logs + "/picard/dedup/{sample}.paired.log"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"]
    shadow: "copy-minimal"
    threads: lambda wildcards, attempt: attempt*2
    resources:
        time=lambda wildcards, attempt: attempt*1440
    wrapper:
        "0.84.0/bio/picard/markduplicates"

rule bam_clipoverlap:
    input:
        bam=results+"/mapping/dedup/{sample}.paired.rmdup.bam",
        ref=REF
    output:
        bam=results+"/mapping/dedup/{sample}.clipped.rmdup.bam",
        met=results+"/qc/clipbam/{sample}.clip.metrics"
    log:
        logs+"/fgbio/clipbam/{sample}.log"
    conda:
        # "../envs/bamutil.yaml"
        "../envs/fgbio.yaml"
    shadow: "copy-minimal"
    resources:
        time=lambda wildcards, attempt: attempt*1440
    shell:
        """
        samtools sort -n -u {input.bam} | fgbio ClipBam -i /dev/stdin \
            -o {output.bam} -r {input.ref} -m {output.met} \
            --clip-overlapping-reads=true -S Coordinate 2> {log}
        """

rule dedup_merged:
    input:
        results+"/mapping/mapped/{sample}.merged.bam"
    output:
        json=results+"/mapping/dedup/{sample}.merged.dedup.json",
        hist=results+"/mapping/dedup/{sample}.merged.hist",
        log=results+"/mapping/dedup/{sample}.merged.log",
        bam=temp(results+"/mapping/dedup/{sample}.merged_rmdup.bam"),
        bamfin=temp(results+"/mapping/dedup/{sample}.merged.rmdup.bam")
    log:
        logs + "/dedup/{sample}.merged.log"
    conda:
        "../envs/dedup.yaml"
    shadow: "copy-minimal"
    params:
        outdir=results+"/mapping/dedup"
    resources:
        time=lambda wildcards, attempt: attempt*1440
    shell:
        """
        dedup -i {input} -m -u -o {params.outdir} 2> {log}
        samtools sort -o {output.bamfin} {output.bam} 2>> {log}
        """

def get_dedup_bam(wildcards):
    # Determines if bam should use Picard or DeDup for duplicate removal
    s = wildcards.sample
    if s in samples.index[samples.depth == "low"]:
        return results+"/mapping/dedup/"+s+".merged.rmdup.bam"
    elif s in samples.index[samples.depth == "high"]:
        return results+"/mapping/dedup/"+s+".clipped.rmdup.bam"

rule finalize_rmdup:
    input:
        get_dedup_bam
    output:
        results+"/mapping/dedup/{sample}.rmdup.bam"
    shell:
        """
        cp {input} {output}
        """

rule realignertargetcreator:
    input:
        bam=results+"/mapping/dedup/{sample}.rmdup.bam",
        bai=results+"/mapping/dedup/{sample}.rmdup.bam.bai",
        ref=REF,
        dic=REF+".dict",
        fai=REF+".fai"
    output:
        intervals=results+"/mapping/indelrealign/{sample}.rmdup.intervals"
    log:
        logs+"/gatk/realignertargetcreator/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shadow: "copy-minimal"
    threads: lambda wildcards, attempt: attempt*2
    resources:
        time=lambda wildcards, attempt: attempt*720
    shell:
        """
        gatk3 -T RealignerTargetCreator -nt {threads} -I {input.bam} \
            -R {input.ref} -o {output.intervals} 2> {log}
        """

rule indelrealigner:
    input:
        bam=results+"/mapping/dedup/{sample}.rmdup.bam",
        bai=results+"/mapping/dedup/{sample}.rmdup.bam.bai",
        intervals=results+"/mapping/indelrealign/{sample}.rmdup.intervals",
        ref=REF,
        dic=REF+".dict",
        fai=REF+".fai"
    output:
        realigned=results+"/mapping/{sample}.rmdup.realn.bam"
    log:
        logs+"/gatk/indelrealigner/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shadow: "copy-minimal"
    threads: lambda wildcards, attempt: attempt*4
    resources:
        time=lambda wildcards, attempt: attempt*1440
    shell:
        """
        gatk3 -T IndelRealigner -R {input.ref} -I {input.bam} \
            -targetIntervals {input.intervals} -o {output.realigned} 2> {log}
        """

ruleorder: samtools_index > samtools_index_temp

rule samtools_index_temp:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        results+"/samtools/index_temp/{prefix}.log"
    shadow: "copy-minimal"
    wrapper:
        "0.84.0/bio/samtools/index"

rule samtools_index:
    input:
        results+"/mapping/{sample}{dp}.rmdup.realn.bam"
    output:
        results+"/mapping/{sample}{dp}.rmdup.realn.bam.bai"
    log:
        results+"/samtools/index/{sample}{dp}.log"
    shadow: "copy-minimal"
    wrapper:
        "0.84.0/bio/samtools/index"

rule samtools_subsample:
    input:
        bam=results+"/mapping/{sample}.rmdup.realn.bam"
    output:
        results+"/mapping/{sample}{dp}.rmdup.realn.bam"
    log:
        logs+"/samtools/subsample/{sample}{dp}.log"
    conda:
        "../envs/samtools.yaml"
    shadow: "copy-minimal"
    resources:
        time=lambda wildcards, attempt: attempt*720
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
            ln -sf $(basename {input.bam}) {output}
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
