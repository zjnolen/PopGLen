import pandas as pd
from snakemake.utils import validate, min_version

# Set minimum snakemake version
min_version("6.2.1")

# load config and sample sheets
configfile: "config.yaml"

samples_df = pd.read_table(config["samples"]).set_index("sample", drop=False)
sample_list = list(samples_df['sample'])
pop_list = list(set(samples_df['population']))

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        "data/beagle/2020modern.beagle.gz",
        expand("data/bams/{population}.bamlist", population=pop_list)

rule download_index_ref:
    """
    Downloads reference genome and indexes with samtools and bwa
    """
    output:
        protected("data/reference/20200120.hicanu.purge.prim.fasta.gz"),
        protected("data/reference/20200120.hicanu.purge.prim.fasta.gz.amb"),
        protected("data/reference/20200120.hicanu.purge.prim.fasta.gz.ann"),
        protected("data/reference/20200120.hicanu.purge.prim.fasta.gz.bwt"),
        protected("data/reference/20200120.hicanu.purge.prim.fasta.gz.pac"),
        protected("data/reference/20200120.hicanu.purge.prim.fasta.gz.sa"),
        protected("data/reference/20200120.hicanu.purge.prim.fasta.gz.fai"),
        protected("data/reference/20200120.hicanu.purge.prim.fasta.gz.gzi")
    log:
        "results/logs/prep_reference.log"
    shell:
        """
        mkdir -p data/reference
        cd data/reference

        wget https://darwin.cog.sanger.ac.uk/insects/Polyommatus_icarus/ilPolIcar1/assemblies/working/20200120.hicanu.purge/20200120.hicanu.purge.prim.fasta.gz

        bwa index 20200120.hicanu.purge.prim.fasta.gz
        samtools faidx 20200120.hicanu.purge.prim.fasta.gz
        """

rule adapterremoval:
    """
    Remove adapters and trim low quality bases at the ends of reads
    """
    output:
        protected("data/fastq_adaptrem/{sample_id}.pair1.truncated.gz"),
        protected("data/fastq_adaptrem/{sample_id}.pair2.truncated.gz"),
        protected("data/fastq_adaptrem/{sample_id}.settings"),
        protected("data/fastq_adaptrem/{sample_id}.discarded.gz"),
        protected("data/fastq_adaptrem/{sample_id}.singleton.truncated.gz")
    params:
        ngi_id = lambda wildcards: samples_df['ngi_id'][wildcards.sample_id]
    log:
        "results/logs/remove_adapters/{sample_id}.log"
    resources:
        runtime = lambda wildcards, attempt: attempt*240
    threads: 4
    shell:
        """
        mkdir -p data/fastq_adaptrem

        AdapterRemoval --file1 data/fastq_raw/{params.ngi_id}*_R1_001.fastq.gz \
        --file2 data/fastq_raw/{params.ngi_id}*_R2_001.fastq.gz --gzip \
        --basename data/fastq_adaptrem/{wildcards.sample_id} --trimns \
        --trimqualities --threads {threads}
        """

rule bwa_map:
    """
    Maps reads to reference using bwa-mem
    """
    input:
        reads=["data/fastq_adaptrem/{sample_id}.pair1.truncated.gz",
            "data/fastq_adaptrem/{sample_id}.pair2.truncated.gz"],
        ref="data/reference/20200120.hicanu.purge.prim.fasta.gz",
        refidx="data/reference/20200120.hicanu.purge.prim.fasta.gz.bwt"
    output:
        "data/bams/{sample_id}.sorted.bam"
    log:
        "results/logs/bwa_map/stdout.bwa_map.{sample_id}.log"
    params:
        extra=r"-R '@RG\tID:2020-01\tSM:{sample_id}\tPL:ILLUMINA'",
        sort="samtools",
        sort_order="coordinates"
    resources:
        runtime = 1440
    threads: 20
    shell:
        """
        (bwa mem -t {threads} {input.ref} {input.reads} | samtools sort \
            -o {output[0]}) > {log}
        """

rule picard_dedup:
    """
    Removes duplicate reads from bam files
    """
    input:
        "data/bams/{sample_id}.sorted.bam"
    output:
        bam=protected("data/bams/{sample_id}.sorted.dedup.bam"),
        metrics=protected("data/bams/{sample_id}.sorted.dedup.metrics.txt")
    log:
        "results/logs/picard_dedup/stdout.picard_dedup.{sample_id}.log"
    resources:
        runtime = 240
    shell:
        """
        (picard MarkDuplicates REMOVE_DUPLICATES=true I={input[0]} O={output.bam} \
            M={output.metrics}) > {log}
        """

rule samtools_index_bam:
    """
    Indexes bam file after sorting and duplicate removal
    """
    input:
        "data/bams/{sample_id}.sorted.dedup.bam"
    output:
        "data/bams/{sample_id}.sorted.dedup.bam.bai"
    log:
        "results/logs/samtools_index_bam/stdout.samtools_index_bam.{sample_id}.log"
    resources:
        runtime = 10
    shell:
        """
        (samtools index {input[0]}) > {log}
        """

rule bamlist_all:
    """
    Make bamlist containing all individuals
    """
    input:
        bams=expand("data/bams/{sample_id}.sorted.dedup.bam",
            sample_id = sample_list),
        bais=expand("data/bams/{sample_id}.sorted.dedup.bam.bai",
            sample_id = sample_list)
    output:
        protected("data/bams/2020modern.bamlist")
    resources:
        runtime = 5
    shell:
        """
        (readlink -f {input.bams} | perl -pe 'chomp if eof') > {output}
        """

rule angsd_beagle_all:
    """
    Make a beagle file from all samples
    """
    input:
        "data/bams/2020modern.bamlist",
        "data/reference/20200120.hicanu.purge.prim.fasta.gz"
    output:
        "data/beagle/2020modern.beagle.gz"
    params:
        outpre="data/beagle/2020modern"
    log:
        "results/logs/angsd_beagle_all/stdout.angsd_beagle_all.log"
    resources:
        runtime = 14400
    threads: 10
    shell:
        """
        (angsd -GL 1 -nThreads {threads} -doGlf 2 -doMajorMinor 1 -minMapQ 30 \
            -c 50 -uniqueOnly 1 -minQ 20 -baq 1 -doMaf 1 -SNP_pval 2e-6 \
            -remove_bads 1 -minInd 36 -bam {input[0]} -out {params.outpre} \
            -ref data/reference/20200120.hicanu.purge.prim.fasta.gz) > \
            {log}
        """

def population2samples(wildcards):
    return expand(
        "data/bams/{sample_id}.sorted.dedup.bam", sample_id =
            samples_df.index[samples_df.population == wildcards.population]
    )

rule bamlist_perpop:
    """
    Make a bamlist from all samples in each population
    """
    input:
        bams = population2samples
    output:
        "data/bams/{population}.bamlist"
    resources:
        runtime = 5
    shell:
        """
        (readlink -f {input.bams} | perl -pe 'chomp if eof') > {output}
        """

rule angsd_perpop_saf:
    """
    Make a saf file from all samples in each population
    """
    input:
        "data/bams/{population}.bamlist",
        "data/reference/20200120.hicanu.purge.prim.fasta.gz"
    output:
        "data/saf/{population}.saf.gz"
    params:
        outpre="data/saf/{population}"
    log:
        "results/logs/angsd_perpop_saf/stdout.angsd_perpop_saf.{population}.log"
    resources:
        runtime = 14400
    threads: 10
    shell:
        """
        (angsd -GL 1 -nThreads {threads} -dosaf 1 -doMajorMinor 1 \
            -c 50 -uniqueOnly 1 -minQ 20 -baq 1 -doMaf 1 -SNP_pval 2e-6 \
            -remove_bads 1 -minInd 7 -bam {input[0]} -out {params.outpre} \
            -ref data/reference/20200120.hicanu.purge.prim.fasta.gz) > {log}
        """
