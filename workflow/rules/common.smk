from snakemake.utils import validate
import pandas as pd
import os

# load and validate config file
if os.path.exists("config/config.yaml"):
    configfile: "config/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml", set_default=True)

# set results directory
results = config["paths"]["results"]

# load and validate sample sheet
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(df, schema="../schemas/samples.schema.yaml")

# load and validate unit sheet
units = pd.read_table(config["units"]).set_index("sample", drop=False)

##### Helper functions #####

def genome_file():
    if pd.isna(config['reference']['fasta_path']):
        fasta = config['reference']['fasta_url']
    else:
        fasta = config['reference']['fasta_path']
    if fasta.endswith('.gz'):
        return "resources/reference/genome.fa.gz"
    else:
        return "resources/reference/genome.fa"

def get_raw_fastq(wildcards):
    # checks if raw sequencing data is specified locally and sets fastp input
    # as either raw local files or to download from SRA
    unit = units.loc[wildcards.sample, ["fq1", "fq2"]]
    return [unit.fq1,unit.fq2]

def get_trimmed_reads(wildcards):
    return rules.fastp_pe.output.merged

def get_read_group(wildcards):
    return r"-R '@RG\tID:{unit}\tSM:{sample}\tPL:{platform}'".format(
        unit=units.loc[wildcards.sample, "unit"],
        sample=wildcards.sample,
        platform=units.loc[wildcards.sample, "platform"]
    )