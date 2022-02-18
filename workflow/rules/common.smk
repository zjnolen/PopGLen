from snakemake.utils import validate
import pandas as pd
import os
import urllib
import re

# load and validate config file
if os.path.exists("config/config.yaml"):
    configfile: "config/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml", set_default=True)

# set results directory
results = config["paths"]["results"]+"/"+config["dataset"]

# set intermediate directory
intermediate = config["paths"]["intermediate"]+"/"+config["dataset"]

# load and validate sample sheet
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(df, schema="../schemas/samples.schema.yaml")

# load and validate unit sheet
units = pd.read_table(config["units"]).set_index("sample", drop=False)

##### Helper functions #####

######### Mapping ##########

def genome_file():
    # Returns path and name of reference after processing
    fasta = config['reference']['fasta']

    if os.path.isfile(fasta):
        fasta = fasta
    else:
        fasta = "resources/reference/"+os.path.basename(fasta)
    
    if fasta.endswith('.gz'):
        return os.path.splitext(fasta)[0]

def get_raw_fastq(wildcards):
    # Checks if raw sequencing data is specified locally and sets fastp 
    # input as either raw local files or to download from SRA
    unit = units.loc[wildcards.sample, ["fq1", "fq2"]]
    return [unit.fq1,unit.fq2]

def get_read_group(wildcards):
    return r"-R '@RG\tID:{unit}\tSM:{sample}\tLB:{sample}\tPL:{platform}'"\
        .format(
            unit=units.loc[wildcards.sample, "unit"],
            sample=wildcards.sample,
            platform=units.loc[wildcards.sample, "platform"]
        )

def get_subsample_prop(wildcards):
    cov_histo = pd.read_table(results + \
        "/qualimap/{wildcards.sample}_mem_dedup/raw_data_qualimapReport/"\
        "coverage_histogram.txt")

    cov = sum(cov_histo["#Coverage"] \
            * cov_histo["Number of genomic locations"]) \
        / sum(cov_histo["Number of genomic locations"])
    
    return "-s " + str(1 / (cov / float({wildcards.cov})))

########## ANGSD ###########

def get_bamlist_bams(wildcards):
    # Checks if rule is looking for a population or a sample by 
    # comparing wildcard value with population list
    pop = wildcards.population
    if re.search("subcov", pop):
        rawpop = pop.rsplit('_', 1)[0]
        suffix = "_" + pop.split('_')[-1]
    else:
        rawpop = pop
        suffix = ""
    if rawpop == config["dataset"]:
        return expand(results + "/dedup/{sample}" + suffix + ".bam", sample
            = samples.index)
    elif rawpop in samples.population.values and pop not in samples.index:
        return expand(results + "/dedup/{sample}" + suffix + ".bam", sample
            = samples.index[samples.population == wildcards.population])
    elif rawpop in samples.index and pop not in samples.population.values:
        return results + "/dedup/{population}.bam"
    elif rawpop in samples.index and pop in samples.population.values:
        print("ERROR: Please ensure sample names do not match population " \
            "names")

def get_bamlist_bais(wildcards):
    pop = wildcards.population
    if re.search("subcov", pop):
        rawpop = pop.rsplit('_', 1)[0]
        suffix = "_" + pop.split('_')[-1]
    else:
        rawpop = pop
        suffix = ""
    if rawpop == config["dataset"]:
        return expand(results + "/dedup/{sample}" + suffix + ".bam.bai", sample
            = samples.index)
    elif rawpop in samples.population.values and pop not in samples.index:
        return expand(results + "/dedup/{sample}" + suffix + ".bam.bai", sample
            = samples.index[samples.population == wildcards.population])
    elif rawpop in samples.index and pop not in samples.population.values:
        return results + "/dedup/{population}.bam.bai"
    elif rawpop in samples.index and pop in samples.population.values:
        print("ERROR: Please ensure sample names do not match population " \
            "names")

def aggregate_beagles(wildcards):
    checkpoint_output = checkpoints.chromosome_list.get(**wildcards).output[0]
    chroms = pd.read_table(checkpoint_output, header=None)
    excl_chr = [config["reference"]["mito"]] + config["reference"]["excl_chr"]
    chroms = chroms[~chroms[0].isin(excl_chr)]
    return expand(results + "/angsd/beagle/" + wildcards.population + 
                "_{chrom}_md" + wildcards.miss + ".beagle.gz", chrom=chroms[0])

def aggregate_ngsLD(wildcards):
    checkpoint_output = checkpoints.chromosome_list.get(**wildcards).output[0]
    chroms = pd.read_table(checkpoint_output, header=None)
    excl_chr = [config["reference"]["mito"]] + config["reference"]["excl_chr"]
    chroms = chroms[~chroms[0].isin(excl_chr)]
    return expand(results + "/ngsLD/" + wildcards.population + "_{chrom}_md"
                + wildcards.miss + "_pruned.sites", chrom=chroms[0])