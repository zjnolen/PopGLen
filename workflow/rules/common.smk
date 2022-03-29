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

# set logs directory
logs = config["paths"]["logs"]+"/"+config["dataset"]

# load and validate sample sheet
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(df, schema="../schemas/samples.schema.yaml")

# load and validate unit sheet
units = pd.read_table(config["units"]).set_index("sample", drop=False)

# get a list of all the populations
if config["populations"]:
    pop_list = config["populations"]
else:
    pop_list = samples.population.unique().tolist()

##### Helper functions #####

#### Mapping & Reference ###

def genome_file():
    # Returns path and name of reference after processing
    fasta = config['reference']['fasta']

    if os.path.isfile(fasta):
        fasta = fasta
    else:
        fasta = "resources/reference/"+os.path.basename(fasta)
    
    if fasta.endswith('.gz'):
        return os.path.splitext(fasta)[0]
    else:
        return fasta

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

def get_autos():
    #set up variables for contig filtering
    mito_sex = config["reference"]["mito"] + config["reference"]["sex"]
    if config["reference"]["min_size"] > 0:
        min_size = config["reference"]["min_size"]
    else:
        min_size = 0
    
    #get autosomes from method specified in config
    if config["reference"]["autosome_file"]:
        auto_list = pd.read_table(
                        config["reference"]["autosome_file"],
                        header=None, squeeze=True
                        )
        auto_list = auto_list[~auto_list.isin(mito_sex)]
        autosomes = auto_list.tolist()
    elif config["reference"]["autosome_list"]:
        auto_list = pd.Series(config["reference"]["autosome_list"])
        auto_list = auto_list[~auto_list.isin(mito_sex)]
        autosomes = auto_list.tolist()
    else:
        with checkpoints.samtools_faidx.get().output[0].open() as fai:
            contigs = pd.read_table(
                            fai, header = None, usecols=[0,1]
                            )
        if config["reference"]["exclude_file"]:
            excl = pd.read_table(
                        config["reference"]["exclude_file"],
                        header=None, squeeze=True
                        ).tolist()
        elif config["reference"]["exclude_list"]:
            excl = config["reference"]["exclude_list"]
        else:
            excl = []
        excl = excl + mito_sex
        contigs = pd.Series(contigs[contigs[1] > min_size][0])
        contigs = contigs[~contigs.isin(excl)]
        autosomes = contigs.tolist()
    
    return autosomes

########## ANGSD ###########

def get_miss_data_prop(wildcards):
    pop = wildcards.population
    if pop == config["dataset"]:
        return config["params"]["angsd"]["max_missing_dataset"]
    elif pop in samples.population.values:
        return config["params"]["angsd"]["max_missing_pop"]
    elif pop in samples.index:
        return 0.0

def get_samples_from_pop(population):
    pop = population
    if pop == config["dataset"]:
        return samples.index.values.tolist()
    elif pop in samples.population.values and pop not in samples.index:
        return samples.index[samples.population == pop].values.tolist()
    elif pop in samples.index and pop not in samples.population.values:
        return [pop]
    elif pop in samples.index and pop in samples.population.values:
        print("ERROR: Ensure no sample shares a name with a population.")
    else:
        print("ERROR: Population queried does not exist in dataset.")

def get_bamlist_bams(wildcards):
    pop = wildcards.population
    return expand(results + "/dedup/{sample}.bam", 
                sample = get_samples_from_pop(pop))

def get_bamlist_bais(wildcards):
    pop = wildcards.population
    return expand(results + "/dedup/{sample}.bam.bai", 
                sample = get_samples_from_pop(pop))

def get_intersect_inputs(wildcards):
    if wildcards.type == "population":
        return expand(results+"/sites/{group}_chr{chrom}_filtautos.bed", group=pop_list, chrom=wildcards.chrom)
    elif wildcards.type == "sample":
        return expand(results+"/sites/{group}_chr{chrom}_filtautos.bed", group=samples.index, chrom=wildcards.chrom)

def get_sites_file(wildcards):
    if wildcards.sites == ".intersect":
        pop = wildcards.population
        if pop in samples.population.values and pop not in samples.index:
            return results+"/sites/population_autos_intersect.sites.idx"
        elif pop in samples.index and pop not in samples.population.values:
            return results+"/sites/sample_autos_intersect.sites.idx"
        elif pop in samples.index and pop in samples.population.values:
            print("ERROR: Ensure no sample shares a name with a population.")
        else:
            print("ERROR: Population queried does not exist in dataset.")

def get_sites_option(wildcards):
    if wildcards.sites == ".intersect":
        pop = wildcards.population
        if pop in samples.population.values and pop not in samples.index:
            return "-sites "+results+"/sites/population_autos_intersect.sites"
        elif pop in samples.index and pop not in samples.population.values:
            return "-sites "+results+"/sites/sample_autos_intersect.sites"
        elif pop in samples.index and pop in samples.population.values:
            print("ERROR: Ensure no sample shares a name with a population.")
        else:
            print("ERROR: Population queried does not exist in dataset.")