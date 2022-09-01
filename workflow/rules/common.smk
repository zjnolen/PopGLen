from snakemake.utils import validate
import pandas as pd
import os
import itertools
import csv

# Set up software containers (for easier version updating)
angsd_container="docker://zjnolen/angsd:0.938"
pcangsd_container="docker://zjnolen/pcangsd:1.10"

# load and validate config file
# if os.path.exists("config/config.yaml"):
#     configfile: "config/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml", set_default=True)

dataset = config["dataset"]

# set results directory
results = "results/"+dataset

# set intermediate directory
intermediate = "intermediate/"+dataset

# set logs directory
logs = "logs/"+dataset

# set up reference variables
REF_NAME=config["reference"]["name"]
REF_DIR=os.path.dirname(config["reference"]["fasta"])
REF=config["reference"]["fasta"]

# define genome chunks to break up analysis
def chunkify(reference_fasta, chunk_size):
    contigs = []
    with open(reference_fasta, "r") as fasta_in:
        for header, seq in itertools.groupby(fasta_in, lambda x: x.startswith(">")):
            if header:
                contig = next(seq).strip(">").strip()
            seq_length = len("".join(seq).replace("\n",""))
            if seq_length > 0:
                contigs = contigs + [[contig,seq_length]]
    df = pd.DataFrame(contigs, columns = ['contig','length']).set_index('contig')
    df.drop(config["reference"]["mito"]+config["reference"]["XZ"]+config["reference"]["exclude"])
    dfs = []
    included = []
    total = 0
    if config["reference"]["min_size"]:
        minsize = config["reference"]["min_size"]
    else:
        minsize = 0
    for n, (index, row) in enumerate(df.iterrows()):
        
        if row['length'] >= minsize:
            total += row['length']
            if total > chunk_size or n+1 >= len(df):
                new_df = pd.DataFrame(included)
                dfs.append(new_df)
                included = []
                #new_df = df.iloc[0:0, :].copy()
                total = row['length']
            included.append(row)
    
    return dfs

chunks = chunkify(REF, config["chunk_size"])
chunklist = list(range(1,len(chunks)+1))

# load and validate sample sheet
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# drop samples specified in config
samples = samples.drop(config["exclude_ind"])

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

# def all_chroms():
#     with checkpoints.samtools_faidx.get().output[0].open() as fai:
#         chroms = pd.read_table(
#                     fai, 
#                     header=None,
#                     usecols=[0]
#                     ).squeeze("columns")
#     return chroms

# def filt_chroms():
#     with checkpoints.combine_beds.get().output.chr.open() as chroms:
#         return pd.read_table(chroms,header=None,usecols=[0]).squeeze("columns")

# def get_autos():
#     #set up variables for contig filtering
#     mito_sex = config["reference"]["mito"] + config["reference"]["sex"]
#     if config["reference"]["min_size"] > 0:
#         min_size = config["reference"]["min_size"]
#     else:
#         min_size = 0
    
#     #get autosomes from method specified in config
#     if config["reference"]["autosome_file"]:
#         auto_list = pd.read_table(
#                         config["reference"]["autosome_file"],
#                         header=None
#                         ).squeeze("columns")
#         auto_list = auto_list[~auto_list.isin(mito_sex)]
#         autosomes = auto_list.tolist()
#     elif config["reference"]["autosome_list"]:
#         auto_list = pd.Series(config["reference"]["autosome_list"])
#         auto_list = auto_list[~auto_list.isin(mito_sex)]
#         autosomes = auto_list.tolist()
#     else:
#         with checkpoints.samtools_faidx.get().output[0].open() as fai:
#             contigs = pd.read_table(
#                             fai, header = None, usecols=[0,1]
#                             )
#         if config["reference"]["exclude_file"]:
#             excl = pd.read_table(
#                         config["reference"]["exclude_file"],
#                         header=None
#                         ).squeeze("columns").tolist()
#         elif config["reference"]["exclude_list"]:
#             excl = config["reference"]["exclude_list"]
#         else:
#             excl = []
#         excl = excl + mito_sex
#         contigs = pd.Series(contigs[contigs[1] > min_size][0])
#         contigs = contigs[~contigs.isin(excl)]
#         autosomes = contigs.tolist()
    
#     return autosomes

########## ANGSD ###########

def get_miss_data_prop(wildcards):
    pop = wildcards.population
    if pop == "all":
        return config["params"]["angsd"]["max_missing_dataset"]
    elif pop in samples.population.values:
        return config["params"]["angsd"]["max_missing_pop"]
    elif pop in samples.index:
        return 0.0

# The following function is useful for now, but not robust to duplicate 
# names across column types. Needs improvement

def get_samples_from_pop(population):
    pop = population
    if pop == "all":
        return samples.index.values.tolist()
    elif pop in samples.depth.values:
        return samples.index[samples.depth == pop].values.tolist()
    elif pop in samples.population.values and pop not in samples.index:
        return samples.index[samples.population == pop].values.tolist()
    elif pop in samples.index and pop not in samples.population.values:
        return [pop]

def get_bamlist_bams(wildcards):
    pop = wildcards.population
    return expand(results+"/mapping/{sample}{{dp}}.rmdup.realn.bam", 
                sample = get_samples_from_pop(pop))

def get_bamlist_bais(wildcards):
    pop = wildcards.population
    return expand(results+"/mapping/{sample}{{dp}}.rmdup.realn.bam.bai", 
                sample = get_samples_from_pop(pop))

# def get_intersect_inputs(wildcards):
#     if wildcards.type == "population":
#         return expand(results +
#                       "/sites/{group}{{dp}}_chr{{chrom}}_filtautos.bed",
#                       group=pop_list
#                       )
#     elif wildcards.type == "sample":
#         return expand(results +
#                       "/sites/{group}{{dp}}_chr{{chrom}}_filtautos.bed",
#                       group=samples.index
#                       )

# def get_sites_file(wildcards):
#     if wildcards.sites == ".intersect":
#         pop = wildcards.population
#         if pop in samples.population.values and pop not in samples.index:
#             return results+"/sites/population"+wildcards.dp+ \
#                         "_autos_intersect.sites.idx"
#         elif pop in samples.index and pop not in samples.population.values:
#             return results+"/sites/sample"+wildcards.dp+ \
#                         "_autos_intersect.sites.idx"
#         elif pop in samples.index and pop in samples.population.values:
#             print("ERROR: Ensure no sample shares a name with a population.")
#         else:
#             print("ERROR: Population queried does not exist in dataset.")
#     return genome_file()

# def get_sites_option(wildcards):
#     if wildcards.sites == ".intersect":
#         pop = wildcards.population
#         if pop in samples.population.values and pop not in samples.index:
#             return "-sites "+results+"/sites/population"+wildcards.dp+ \
#                         "_autos_intersect.sites"
#         elif pop in samples.index and pop not in samples.population.values:
#             return "-sites "+results+"/sites/sample"+wildcards.dp+ \
#                         "_autos_intersect.sites"
#         elif pop in samples.index and pop in samples.population.values:
#             print("ERROR: Ensure no sample shares a name with a population.")
#         else:
#             print("ERROR: Population queried does not exist in dataset.")
#     else:
#         return ""