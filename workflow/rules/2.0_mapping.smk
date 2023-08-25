# Rules for mapping processed reads to the reference and processing bam files


rule bwa_mem_merged:
    """Map collapsed read pairs for historical samples to reference genome"""
    input:
        reads="results/preprocessing/fastp/{sample}.merged.fastq.gz",
        ref="results/ref/{ref}/{ref}.fa",
        idx=rules.bwa_index.output,
    output:
        temp("results/mapping/mapped/{sample}.{ref}.merged.bam"),
    log:
        "logs/mapping/bwa_mem/{sample}.{ref}.merged.log",
    benchmark:
        "benchmarks/mapping/bwa_mem/{sample}.{ref}.merged.log"
    params:
        extra=lambda w: f"{get_read_group(w)}",
        sort="samtools",
    threads: lambda wildcards, attempt: attempt * 10
    resources:
        runtime=lambda wildcards, attempt: attempt * 2880,
    wrapper:
        "v2.6.0/bio/bwa-mem2/mem"


rule bwa_mem_paired:
    """Map trimmed paired reads from modern samples to reference genome"""
    input:
        reads=expand(
            "results/preprocessing/fastp/{{sample}}.{read}.fastq.gz", read=["R1", "R2"]
        ),
        ref="results/ref/{ref}/{ref}.fa",
        idx=rules.bwa_index.output,
    output:
        bam=temp("results/mapping/mapped/{sample}.{ref}.paired.bam"),
    log:
        "logs/mapping/bwa_mem/{sample}.{ref}.paired.log",
    benchmark:
        "benchmarks/mapping/bwa_mem/{sample}.{ref}.paired.log"
    params:
        extra=lambda w: f"{get_read_group(w)}",
        sort="samtools",
    threads: lambda wildcards, attempt: attempt * 10
    resources:
        runtime=lambda wildcards, attempt: attempt * 2880,
    wrapper:
        "v2.6.0/bio/bwa-mem2/mem"


rule mark_duplicates:
    """Remove duplicate reads from paired end bam files"""
    input:
        bams="results/mapping/mapped/{sample}.{ref}.paired.bam",
    output:
        bam=temp("results/mapping/dedup/{sample}.{ref}.paired.rmdup.bam"),
        metrics="results/mapping/qc/mark_duplicates/{sample}.{ref}.picard.metrics",
    log:
        "logs/mapping/picard/dedup/{sample}.{ref}.paired.log",
    benchmark:
        "benchmarks/mapping/picard/dedup/{sample}.{ref}.paired.log"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"],
    shadow:
        "minimal"
    threads: lambda wildcards, attempt: attempt * 4
    resources:
        # can be memory intensive for big bam files, look into ways of 
        # allocating memory that will work on multiple clusters
        runtime=lambda wildcards, attempt: attempt * 1440,
    wrapper:
        "v1.17.2/bio/picard/markduplicates"


rule bam_clipoverlap:
    """Clip overlapping reads in paired end bam files"""
    input:
        bam="results/mapping/dedup/{sample}.{ref}.paired.rmdup.bam",
        ref="results/ref/{ref}/{ref}.fa",
    output:
        bam=temp("results/mapping/dedup/{sample}.{ref}.clipped.rmdup.bam"),
        bai=temp("results/mapping/dedup/{sample}.{ref}.clipped.rmdup.bam.bai"),
        met="results/mapping/qc/fgbio_clipbam/{sample}.{ref}.fgbio_clip.metrics",
    log:
        "logs/mapping/fgbio/clipbam/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/fgbio/clipbam/{sample}.{ref}.log"
    conda:
        "../envs/fgbio.yaml"
    shadow:
        "minimal"
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 1440,
    shell:
        """
        (samtools sort -n -T {resources.tmpdir} -o {input.bam}.namesort.bam {input.bam}
        fgbio -Xmx{resources.mem_mb}m ClipBam \
                -i {input.bam}.namesort.bam \
                -r {input.ref} -m {output.met} \
                --clip-overlapping-reads=true \
                -o {output.bam}.namesort.bam
        samtools sort -T {resources.tmpdir} -o {output.bam} {output.bam}.namesort.bam
        samtools index {output.bam}) 2> {log}
        """


rule dedup_merged:
    """Remove duplicates from collapsed read bam files"""
    input:
        "results/mapping/mapped/{sample}.{ref}.merged.bam",
    output:
        json="results/mapping/qc/dedup/{sample}.{ref}.dedup.json",
        hist="results/mapping/qc/dedup/{sample}.{ref}.dedup.hist",
        log="results/mapping/qc/dedup/{sample}.{ref}.dedup.log",
        bam=temp("results/mapping/dedup/{sample}.{ref}.merged_rmdup.bam"),
        bamfin=temp("results/mapping/dedup/{sample}.{ref}.merged.rmdup.bam"),
        bai=temp("results/mapping/dedup/{sample}.{ref}.merged.rmdup.bam.bai"),
    log:
        "logs/mapping/dedup/{sample}.{ref}.merged.log",
    benchmark:
        "benchmarks/mapping/dedup/{sample}.{ref}.merged.log"
    conda:
        "../envs/dedup.yaml"
    shadow:
        "minimal"
    params:
        outdir=lambda w, output: os.path.dirname(output.bamfin),
    resources:
        runtime=lambda wildcards, attempt: attempt * 1440,
    shell:
        """
        (dedup -i {input} -m -u -o {params.outdir}
        samtools sort -T {resources.tmpdir} -o {output.bamfin} {output.bam}
        samtools index {output.bamfin}
        mv {params.outdir}/{wildcards.sample}.{wildcards.ref}.merged.dedup.json \
            {output.json}
        mv {params.outdir}/{wildcards.sample}.{wildcards.ref}.merged.hist \
            {output.hist}
        mv {params.outdir}/{wildcards.sample}.{wildcards.ref}.merged.log \
            {output.log}) &> {log}
        """


rule realignertargetcreator:
    """Create intervals database for GATK Indel Realigner"""
    input:
        bam=get_dedup_bam,
        ref="results/ref/{ref}/{ref}.fa",
        dic="results/ref/{ref}/{ref}.dict",
        fai="results/ref/{ref}/{ref}.fa.fai",
    output:
        intervals="results/mapping/indelrealign/{sample}.{ref}.rmdup.intervals",
    log:
        "logs/mapping/gatk/realignertargetcreator/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/gatk/realignertargetcreator/{sample}.{ref}.log"
    conda:
        "../envs/gatk.yaml"
    shadow:
        "minimal"
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        gatk3 -T RealignerTargetCreator -nt {threads} -I {input.bam[0]} \
            -R {input.ref} -o {output.intervals} 2> {log}
        """


rule indelrealigner:
    """Realign reads around indels"""
    input:
        bam=get_dedup_bam,
        intervals="results/mapping/indelrealign/{sample}.{ref}.rmdup.intervals",
        ref="results/ref/{ref}/{ref}.fa",
        dic="results/ref/{ref}/{ref}.dict",
        fai="results/ref/{ref}/{ref}.fa.fai",
    output:
        realigned="results/mapping/bams/{sample}.{ref}.rmdup.realn.bam",
    log:
        "logs/mapping/gatk/indelrealigner/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/gatk/indelrealigner/{sample}.{ref}.log"
    conda:
        "../envs/gatk.yaml"
    shadow:
        "minimal"
    threads: lambda wildcards, attempt: attempt * 4
    resources:
        runtime=lambda wildcards, attempt: attempt * 1440,
    shell:
        """
        gatk3 -Xmx{resources.mem_mb}m -T IndelRealigner -R {input.ref} \
            -I {input.bam[0]} -targetIntervals {input.intervals} \
            -o {output.realigned} 2> {log}
        """


rule samtools_index:
    """Index bam files """
    input:
        "results/mapping/bams/{prefix}.bam",
    output:
        "results/mapping/bams/{prefix}.bam.bai",
    log:
        "logs/mapping/samtools/index/{prefix}.log",
    benchmark:
        "benchmarks/mapping/samtools/index/{prefix}.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input} 2> {log}
        """


rule symlink_bams:
    """
    Link bam files to be used in a dataset into the dataset folder. Leaving bam files 
    outside the dataset folder allows them to be reused in other datasets if processed 
    the same way
    """
    input:
        unpack(get_final_bam),
    output:
        bam="results/datasets/{dataset}/bams/{sample}.{ref}.bam",
        bai="results/datasets/{dataset}/bams/{sample}.{ref}.bam.bai",
    log:
        "logs/{dataset}/symlink_bams/{sample}.{ref}.log",
    conda:
        "../envs/shell.yaml"
    shell:
        """
        ln -sr {input.bam} {output.bam}
        ln -sr {input.bai} {output.bai}
        """


rule samtools_subsample:
    """
    Subsample all bam files down to the same coverage to examine effects of variance in 
    coverage
    """
    input:
        bam=get_bamlist_bams,
        bai=get_bamlist_bais,
        depth="results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}_{population}.depth.sum",
        bed="results/datasets/{dataset}/genotyping/filters/beds/{dataset}_filts.bed",
    output:
        bam="results/datasets/{dataset}/bams/{population}{dp}.bam",
        bai="results/datasets/{dataset}/bams/{population}{dp}.bam.bai",
    log:
        "logs/mapping/samtools/subsample/{dataset}_{population}{dp}.log",
    benchmark:
        "benchmarks/mapping/samtools/subsample/{dataset}_{population}{dp}.log"
    conda:
        "../envs/samtools.yaml"
    shadow:
        "minimal"
    params:
        dp=config["downsample_cov"],
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        dp=$(awk '{{print $2}}' {input.depth})

        prop=$(echo "$dp {params.dp}" | awk '{{print 1 / ($1 / $2)}}')

        if [ `awk 'BEGIN {{print ('$prop' <= 1.0)}}'` = 1 ]; then
            propdec=$(echo $prop | awk -F "." '{{print $2}}')
            samtools view -h -s ${{RANDOM}}.${{propdec}} -q 30 -L {input.bed} -@ {threads} \
                -b {input.bam} > {output.bam} 2> {log}
            samtools index {output.bam} 2>> {log}
        else
            echo "WARNING: Depth of sample is lower than subsample depth." \
                &> {log}
            echo "Subsampled bam will be symlink of original." &>> {log}
            ln -sf $(basename {input.bam}) {output.bam} 2>> {log}
            ln -sf $(basename {input.bai}) {output.bai} 2>> {log}
        fi
        """
