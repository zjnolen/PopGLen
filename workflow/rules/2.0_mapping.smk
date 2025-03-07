# Rules for mapping processed reads to the reference and processing bam files


localrules:
    symlink_bams,
    symlink_userbams,


rule bwa_aln_merged:
    """
    Align reads to reference using BWA ALN algorithm.
    """
    input:
        fastq="results/preprocessing/fastp/{sample}_{unit}_{lib}.merged.fastq.gz",
        idx=rules.bwa_index.output,
    output:
        temp("results/mapping/mapped/{sample}_{unit}_{lib}.{ref}.aln.merged.sai"),
    log:
        "logs/mapping/bwa_aln/{sample}_{unit}_{lib}.{ref}.merged.log",
    benchmark:
        "benchmarks/mapping/bwa_aln/{sample}_{unit}_{lib}.{ref}.merged.log"
    params:
        extra=config["params"]["bwa_aln"]["extra"],
    threads: 20
    resources:
        runtime="7d",
    wrapper:
        "v4.0.0/bio/bwa/aln"


rule bwa_samse_merged:
    """
    Convert alignment into BAM format using bwa samse.
    """
    input:
        fastq="results/preprocessing/fastp/{sample}_{unit}_{lib}.merged.fastq.gz",
        sai="results/mapping/mapped/{sample}_{unit}_{lib}.{ref}.aln.merged.sai",
        idx=rules.bwa_index.output,
    output:
        temp("results/mapping/mapped/{sample}_{unit}_{lib}.{ref}.aln.merged.bam"),
    log:
        "logs/mapping/bwa_samse/{sample}_{unit}_{lib}.{ref}.merged.log",
    benchmark:
        "benchmarks/mapping/bwa_samse/{sample}_{unit}_{lib}.{ref}.merged.log"
    params:
        extra=lambda w: f"-r {get_read_group(w)}",
        sort="samtools",
        sort_order="coordinate",
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime="6h",
    wrapper:
        "v4.0.0/bio/bwa/samse"


rule bwa_mem_paired:
    """Map trimmed paired reads from modern samples to reference genome"""
    input:
        reads=expand(
            "results/preprocessing/fastp/{{sample}}_{{unit}}_{{lib}}.{read}.{{pairing}}.fastq.gz",
            read=["R1", "R2"],
        ),
        ref="results/ref/{ref}/{ref}.fa",
        idx=rules.bwa_index.output,
    output:
        bam=temp("results/mapping/mapped/{sample}_{unit}_{lib}.{ref}.mem.{pairing}.bam"),
    log:
        "logs/mapping/bwa_mem/{sample}_{unit}_{lib}.{ref}.{pairing}.log",
    benchmark:
        "benchmarks/mapping/bwa_mem/{sample}_{unit}_{lib}.{ref}.{pairing}.log"
    wildcard_constraints:
        pairing="paired|uncollapsed",
    params:
        extra=lambda w: f"-R {get_read_group(w)}",
        sorting="samtools",
    threads: lambda wildcards, attempt: attempt * 10
    resources:
        runtime=lambda wildcards, attempt: attempt * 2880,
    wrapper:
        "v4.0.0/bio/bwa/mem"


rule bwa_mem_merged:
    """Map collapsed reads from historical samples to reference genome"""
    input:
        reads="results/preprocessing/fastp/{sample}_{unit}_{lib}.merged.fastq.gz",
        ref="results/ref/{ref}/{ref}.fa",
        idx=rules.bwa_index.output,
    output:
        bam=temp("results/mapping/mapped/{sample}_{unit}_{lib}.{ref}.mem.merged.bam"),
    log:
        "logs/mapping/bwa_mem/{sample}_{unit}_{lib}.{ref}.merged.log",
    benchmark:
        "benchmarks/mapping/bwa_mem/{sample}_{unit}_{lib}.{ref}.merged.log"
    params:
        extra=lambda w: f"-R {get_read_group(w)}",
        sorting="samtools",
    threads: lambda wildcards, attempt: attempt * 10
    resources:
        runtime=lambda wildcards, attempt: attempt * 2880,
    wrapper:
        "v4.0.0/bio/bwa/mem"


rule samtools_merge_collapsed_libs:
    input:
        get_lib_bams,
    output:
        bam=temp("results/mapping/mapped/{sample}_{lib}.{ref}.merged.bam"),
    log:
        "logs/mapping/samtools/merge/{sample}_{lib}.{ref}.merged.log",
    benchmark:
        "benchmarks/mapping/samtools/merge/{sample}_{lib}.{ref}.merged.log"
    threads: lambda wildcards, attempt: attempt * 4
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    wrapper:
        "v4.0.0/bio/samtools/merge"


rule samtools_merge_paired_units:
    input:
        get_paired_bams,
    output:
        bam=temp("results/mapping/mapped/{sample}.{ref}.paired.bam"),
    log:
        "logs/mapping/samtools/merge/{sample}.{ref}.paired.log",
    benchmark:
        "benchmarks/mapping/samtools/{sample}.{ref}.paired.log"
    threads: lambda wildcards, attempt: attempt * 4
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    wrapper:
        "v4.0.0/bio/samtools/merge"


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
        runtime=lambda wildcards, attempt: attempt * 1440,
    wrapper:
        "v4.0.0/bio/picard/markduplicates"


rule dedup_merged:
    """Remove duplicates from collapsed read bam files"""
    input:
        "results/mapping/mapped/{sample}_{lib}.{ref}.merged.bam",
    output:
        json="results/mapping/qc/dedup/{sample}_{lib}.{ref}.dedup.json",
        hist="results/mapping/qc/dedup/{sample}_{lib}.{ref}.dedup.hist",
        log="results/mapping/qc/dedup/{sample}_{lib}.{ref}.dedup.log",
        bam=temp("results/mapping/dedup/{sample}_{lib}.{ref}.merged_rmdup.bam"),
        bamfin=temp("results/mapping/dedup/{sample}_{lib}.{ref}.merged.rmdup.bam"),
        bai=temp("results/mapping/dedup/{sample}_{lib}.{ref}.merged.rmdup.bam.bai"),
    log:
        "logs/mapping/dedup/{sample}_{lib}.{ref}.merged.log",
    benchmark:
        "benchmarks/mapping/dedup/{sample}_{lib}.{ref}.merged.log"
    conda:
        "../envs/dedup.yaml"
    shadow:
        "minimal"
    threads: lambda wildcards, attempt: attempt * 2
    params:
        outdir=lambda w, output: os.path.dirname(output.bamfin),
    resources:
        runtime=lambda wildcards, attempt: attempt * 1440,
    shell:
        """
        (dedup -i {input} -m -u -o {params.outdir}
        samtools sort -T {resources.tmpdir} -o {output.bamfin} {output.bam}
        samtools index {output.bamfin}
        mv {params.outdir}/{wildcards.sample}_{wildcards.lib}.{wildcards.ref}.merged.dedup.json \
            {output.json}
        mv {params.outdir}/{wildcards.sample}_{wildcards.lib}.{wildcards.ref}.merged.hist \
            {output.hist}
        mv {params.outdir}/{wildcards.sample}_{wildcards.lib}.{wildcards.ref}.merged.log \
            {output.log}) &> {log}
        """


rule samtools_merge_dedup:
    """
    Merge deduplicated BAMs of historical sample libraries
    """
    input:
        get_dedup_bams,
    output:
        bam=temp("results/mapping/dedup/{sample}.{ref}.rmdup.bam"),
    log:
        "logs/mapping/samtools/merge/{sample}.{ref}.rmdup.log",
    benchmark:
        "benchmarks/mapping/samtools/{sample}.{ref}.rmdup.log"
    threads: lambda wildcards, attempt: attempt * 4
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    wrapper:
        "v4.0.0/bio/samtools/merge"


rule realignertargetcreator:
    """Create intervals database for GATK Indel Realigner"""
    input:
        bam="results/mapping/dedup/{sample}.{ref}.rmdup.bam",
        bai="results/mapping/dedup/{sample}.{ref}.rmdup.bam.bai",
        ref="results/ref/{ref}/{ref}.fa",
        dict="results/ref/{ref}/{ref}.dict",
        fai="results/ref/{ref}/{ref}.fa.fai",
    output:
        intervals="results/mapping/indelrealign/{sample}.{ref}.rmdup.intervals",
    log:
        "logs/mapping/gatk/realignertargetcreator/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/gatk/realignertargetcreator/{sample}.{ref}.log"
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    wrapper:
        "v4.0.0/bio/gatk3/realignertargetcreator"


rule indelrealigner:
    """Realign reads around indels"""
    input:
        bam="results/mapping/dedup/{sample}.{ref}.rmdup.bam",
        bai="results/mapping/dedup/{sample}.{ref}.rmdup.bam.bai",
        target_intervals="results/mapping/indelrealign/{sample}.{ref}.rmdup.intervals",
        ref="results/ref/{ref}/{ref}.fa",
        dict="results/ref/{ref}/{ref}.dict",
        fai="results/ref/{ref}/{ref}.fa.fai",
    output:
        bam=temp("results/mapping/bams/{sample}.{ref}.rmdup.realn.bam"),
    log:
        "logs/mapping/gatk/indelrealigner/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/gatk/indelrealigner/{sample}.{ref}.log"
    threads: lambda wildcards, attempt: attempt * 4
    resources:
        runtime=lambda wildcards, attempt: attempt * 1440,
    wrapper:
        "v4.0.0/bio/gatk3/indelrealigner"


rule bam_clipoverlap:
    """Clip overlapping reads in paired end bam files"""
    input:
        bam="results/mapping/bams/{sample}.{ref}.rmdup.realn.bam",
        ref="results/ref/{ref}/{ref}.fa",
    output:
        bam="results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.bam",
        log="results/mapping/qc/bamutil_clipoverlap/{sample}.{ref}.rmdup.realn.clipoverlap.stats",
    log:
        "logs/mapping/bamutil/clipoverlap/{sample}.{ref}.rmdup.realn.log",
    benchmark:
        "benchmarks/mapping/bamutil/clipoverlap/{sample}.{ref}.rmdup.realn.log"
    container:
        bamutil_container
    shadow:
        "minimal"
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        bam clipOverlap --in {input.bam} --out {output.bam} --stats 2> {log}
        cat {log} > {output.log}
        """


rule symlink_userbams:
    input:
        bam=lambda w: units.loc[units["sample"] == w.sample, "bam"].values[0],
    output:
        bam="results/mapping/user-provided-bams/{sample}.{ref}.user-processed.bam",
    log:
        "logs/symlink_bams/{sample}.{ref}.user-processed.log",
    container:
        shell_container
    resources:
        runtime="5m",
    shell:
        """
        ln -sr {input.bam} {output.bam}
        """


rule bam_clipoverlap_userbams:
    """Clip overlapping reads in paired end bam files provided by users"""
    input:
        bam="results/mapping/user-provided-bams/{sample}.{ref}.user-processed.bam",
        ref="results/ref/{ref}/{ref}.fa",
    output:
        bam="results/mapping/user-provided-bams/{sample}.{ref}.clip.bam",
        log="results/mapping/user-provided-bams/{sample}.{ref}.clipoverlap.stats",
    log:
        "logs/mapping/bamutil/clipoverlap/{sample}.{ref}.user-processed.log",
    benchmark:
        "benchmarks/mapping/bamutil/clipoverlap/{sample}.{ref}.user-processed.log"
    container:
        bamutil_container
    shadow:
        "minimal"
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        bam clipOverlap --in {input.bam} --out {output.bam} --stats 2> {log}
        cat {log} > {output.log}
        """


ruleorder: symlink_bams > samtools_index
ruleorder: samtools_subsample > samtools_index


rule samtools_index:
    """Index bam files"""
    input:
        "results/{prefix}.bam",
    output:
        "results/{prefix}.bam.bai",
    container:
        samtools_container
    log:
        "logs/mapping/samtools/index/{prefix}.log",
    benchmark:
        "benchmarks/mapping/samtools/index/{prefix}.log"
    resources:
        runtime="1h",
    shell:
        """
        samtools index {input} {output} 2> {log}
        """


rule symlink_bams:
    """
    Link bam files to be used in a dataset into the dataset folder. Leaving bam
    files outside the dataset folder allows them to be reused in other datasets
    if processed the same way
    """
    input:
        unpack(get_final_bam),
    output:
        bam="results/datasets/{dataset}/bams/{sample}.{ref}.bam",
        bai="results/datasets/{dataset}/bams/{sample}.{ref}.bam.bai",
    log:
        "logs/{dataset}/symlink_bams/{sample}.{ref}.log",
    container:
        shell_container
    resources:
        runtime="5m",
    shell:
        """
        ln -sr {input.bam} {output.bam}
        ln -sr {input.bai} {output.bai}
        """


rule samtools_subsample:
    """
    Subsample all bam files down to the same coverage to examine effects of
    variance in coverage
    """
    input:
        bam="results/datasets/{dataset}/bams/{sample}.{ref}.bam",
        bai="results/datasets/{dataset}/bams/{sample}.{ref}.bam.bai",
        depth=depth_file,
    output:
        bam="results/datasets/{dataset}/bams/{sample}.{ref}{dp}.bam",
        bai="results/datasets/{dataset}/bams/{sample}.{ref}{dp}.bam.bai",
    log:
        "logs/mapping/samtools/subsample/{dataset}_{sample}.{ref}{dp}.log",
    benchmark:
        "benchmarks/mapping/samtools/subsample/{dataset}_{sample}.{ref}{dp}.log"
    container:
        samtools_container
    shadow:
        "minimal"
    params:
        dp=lambda w: w.dp.replace(".dp", ""),
        seed=config["params"]["samtools"]["subsampling_seed"],
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        dp=$(awk '{{print $2}}' {input.depth})

        prop=$(echo "$dp {params.dp}" | awk '{{print 1 / ($1 / $2)}}')

        if [ `awk 'BEGIN {{print ('$prop' <= 1.0)}}'` = 1 ]; then
            propdec=$(echo $prop | awk -F "." '{{print $2}}')
            samtools view -h -F 4 -q 30 -@ {threads} -u {input.bam} |
                samtools view -h -s {params.seed}.${{propdec}} -@ {threads} -b \
                > {output.bam} 2> {log}
            samtools index {output.bam} 2>> {log}
        else
            echo "WARNING: Depth of sample is lower than subsample depth." \
                &> {log}
            echo "Subsampled bam will be symlink of original." &>> {log}
            ln -sf $(basename {input.bam}) {output.bam} 2>> {log}
            ln -sf $(basename {input.bai}) {output.bai} 2>> {log}
        fi
        """
