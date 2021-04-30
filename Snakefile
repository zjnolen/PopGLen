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
        expand("data/fastq_adaptrem/{sample_id}.collapsed.gz",
            sample_id = sample_list)

rule remove_adapters:
    """
    Remove adapters and trim low quality bases at the ends of reads
    """
    output:
        "data/fastq_adaptrem/{sample_id}.collapsed.gz"
    params:
        ngi_id = lambda wildcards: samples_df['ngi_id'][wildcards.sample_id]
    log:
        "results/logs/remove_adapters/{sample_id}.log"
    resources:
        runtime_min=60
    threads: 4
    shell:
        """
        mkdir -p data/fastq_adaptrem

        AdapterRemoval --file1 data/fastq_raw/{params.ngi_id}*_R1_001.fastq.gz \
        --file2 data/fastq_raw/{params.ngi_id}*_R2_001.fastq.gz --gzip \
        --basename data/fastq_adaptrem/{wildcards.sample_id} --trimns \
        --trimquals --threads {threads}
        """
