rule bwa_mem_merged:
    input:
        merged=rules.fastp_mergedout.output.merged,
        idx=rules.bwa_index.output
    output:
        bam=temp("results/mapping/mapped/{sample}.merged.bam")
    log:
        "logs/mapping/bwa_mem/{sample}.merged.log"
    params:
        index=REF,
        rg=get_read_group
    conda:
        "../envs/mapping.yaml"
    shadow: "copy-minimal"
    threads: lambda wildcards, attempt: attempt*10
    resources:
        time=lambda wildcards, attempt: attempt*2880
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
        bam=temp("results/mapping/mapped/{sample}.paired.bam")
    log:
        "logs/mapping/bwa_mem/{sample}.paired.log"
    params:
        index=REF,
        rg=get_read_group
    conda:
        "../envs/mapping.yaml"
    shadow: "copy-minimal"
    threads: lambda wildcards, attempt: attempt*10
    resources:
        time=lambda wildcards, attempt: attempt*2880
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
        "results/mapping/mapped/{sample}.paired.bam"
    output:
        bam=temp("results/mapping/dedup/{sample}.paired.rmdup.bam"),
        metrics="results/mapping/qc/mark_duplicates/{sample}.picard.metrics"
    log:
        "logs/mapping/picard/dedup/{sample}.paired.log"
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
        bam="results/mapping/dedup/{sample}.paired.rmdup.bam",
        ref=REF
    output:
        bam=temp("results/mapping/dedup/{sample}.clipped.rmdup.bam"),
        met="results/mapping/qc/fgbio_clipbam/{sample}.fgbio_clip.metrics"
    log:
        "logs/mapping/fgbio/clipbam/{sample}.log"
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
        "results/mapping/mapped/{sample}.merged.bam"
    output:
        json="results/mapping/qc/dedup/{sample}.dedup.json",
        hist="results/mapping/qc/dedup/{sample}.dedup.hist",
        log="results/mapping/qc/dedup/{sample}.dedup.log",
        bam=temp("results/mapping/dedup/{sample}.merged_rmdup.bam"),
        bamfin=temp("results/mapping/dedup/{sample}.merged.rmdup.bam")
    log:
        "logs/mapping/dedup/{sample}.merged.log"
    conda:
        "../envs/dedup.yaml"
    shadow: "copy-minimal"
    params:
        outdir="results/mapping/dedup/"
    resources:
        time=lambda wildcards, attempt: attempt*1440
    shell:
        """
        dedup -i {input} -m -u -o {params.outdir} 2> {log}
        samtools sort -o {output.bamfin} {output.bam} 2>> {log}
        mv {params.outdir}{wildcards.sample}.merged.dedup.json \
            {output.json} 2>> {log}
        mv {params.outdir}{wildcards.sample}.merged.hist \
            {output.hist} 2>> {log}
        mv {params.outdir}{wildcards.sample}.merged.log \
            {output.log} 2>> {log}
        """

def get_dedup_bam(wildcards):
    # Determines if bam should use Picard or DeDup for duplicate removal
    s = wildcards.sample
    if s in samples.index[samples.depth == "low"]:
        return ["results/mapping/dedup/"+s+".merged.rmdup.bam",
                "results/mapping/dedup/"+s+".merged.rmdup.bam.bai"]
    elif s in samples.index[samples.depth == "high"]:
        return ["results/mapping/dedup/"+s+".clipped.rmdup.bam",
                "results/mapping/dedup/"+s+".clipped.rmdup.bam.bai"]

ruleorder: finalize_rmdup > samtools_index_temp

rule finalize_rmdup:
    input:
        get_dedup_bam
    output:
        temp("results/mapping/dedup/{sample}.rmdup.bam"),
        temp("results/mapping/dedup/{sample}.rmdup.bam.bai")
    shell:
        """
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        """

rule realignertargetcreator:
    input:
        bam="results/mapping/dedup/{sample}.rmdup.bam",
        bai="results/mapping/dedup/{sample}.rmdup.bam.bai",
        ref=REF,
        dic=REF+".dict",
        fai=REF+".fai"
    output:
        intervals="results/mapping/indelrealign/{sample}.rmdup.intervals"
    log:
        "logs/mapping/gatk/realignertargetcreator/{sample}.log"
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
        bam="results/mapping/dedup/{sample}.rmdup.bam",
        bai="results/mapping/dedup/{sample}.rmdup.bam.bai",
        intervals="results/mapping/indelrealign/{sample}.rmdup.intervals",
        ref=REF,
        dic=REF+".dict",
        fai=REF+".fai"
    output:
        realigned="results/mapping/bams/{sample}.rmdup.realn.bam"
    log:
        "logs/mapping/gatk/indelrealigner/{sample}.log"
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
        temp("{prefix}.bam.bai")
    log:
        "logs/mapping/samtools/index_temp/{prefix}.log"
    shadow: "copy-minimal"
    wrapper:
        "0.84.0/bio/samtools/index"

rule samtools_index:
    input:
        "results/mapping/bams/{sample}{dp}.rmdup.realn.bam"
    output:
        "results/mapping/bams/{sample}{dp}.rmdup.realn.bam.bai"
    log:
        "logs/mapping/samtools/index/{sample}{dp}.log"
    shadow: "copy-minimal"
    wrapper:
        "0.84.0/bio/samtools/index"

rule samtools_subsample:
    input:
        bam="results/mapping/bams/{sample}.rmdup.realn.bam"
    output:
        "results/mapping/bams/{sample}{dp}.rmdup.realn.bam"
    log:
        "logs/mapping/samtools/subsample/{sample}{dp}.log"
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