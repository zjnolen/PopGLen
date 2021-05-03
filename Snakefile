import pandas as pd
from snakemake.utils import validate, min_version

# Set minimum snakemake version
min_version("6.2.1")

# load config and sample sheets
configfile: "config.yaml"

samples_df = pd.read_table(config["samples"]).set_index("sample", drop=False)
sample_list = list(samples_df['sample'])

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        expand("data/fastq_adaptrem/{sample_id}.pair1.truncated.gz",
            sample_id = sample_list),

rule prep_reference:
    """
    Downloads reference genome and indexes with samtools and bwa
    """
    output:
        "data/reference/20200120.hicanu.purge.prim.fasta.gz",
        "data/reference/20200120.hicanu.purge.prim.fasta.gz.amb",
        "data/reference/20200120.hicanu.purge.prim.fasta.gz.ann",
        "data/reference/20200120.hicanu.purge.prim.fasta.gz.bwt",
        "data/reference/20200120.hicanu.purge.prim.fasta.gz.pac",
        "data/reference/20200120.hicanu.purge.prim.fasta.gz.sa",
        "data/reference/20200120.hicanu.purge.prim.fasta.gz.fai",
        "data/reference/20200120.hicanu.purge.prim.fasta.gz.gzi"
    log:
        "results/logs/prep_reference.log"
    resources:
        runtime = 60
    threads: 1
    shell:
        """
        mkdir -p data/reference
        cd data/reference

        wget https://darwin.cog.sanger.ac.uk/insects/Polyommatus_icarus/ilPolIcar1/assemblies/working/20200120.hicanu.purge/20200120.hicanu.purge.prim.fasta.gz

        bwa index 20200120.hicanu.purge.prim.fasta.gz
        samtools faidx 20200120.hicanu.purge.prim.fasta.gz
        """

rule remove_adapters:
    """
    Remove adapters and trim low quality bases at the ends of reads
    """
    output:
        protected("data/fastq_adaptrem/{sample_id}.pair1.truncated.gz"),
        protected("data/fastq_adaptrem/{sample_id}.pair2.truncated.gz"),
        protected("data/fastq_adaptrem/{sample_id}.settings),
        protected("data/fastq_adaptrem/{sample_id}.discarded.gz),
        protected("data/fastq_adaptrem/{sample_id}.singleton.truncated.gz)
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

