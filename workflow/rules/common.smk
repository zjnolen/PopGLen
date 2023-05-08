from snakemake.utils import validate
import pandas as pd
import os
import itertools
import csv

# Set up software containers (for easier version updating)
angsd_container="docker://zjnolen/angsd:0.938"
pcangsd_container="docker://zjnolen/pcangsd:1.10"
evaladmix_container="docker://zjnolen/evaladmix:0.961"
ngsf_hmm_container="docker://zjnolen/ngsf-hmm:20200722-2df9690"
mapdamage_container="docker://quay.io/biocontainers/mapdamage2:2.2.1--pyr40_0"
ngsrelate_container="docker://zjnolen/ngsrelate:20220925-ec95c8f"

# load and validate config file
# if os.path.exists("config/config.yaml"):
#     configfile: "config/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml", set_default=True)

dataset = config["dataset"]

# set results directory
results = "results/results/"+dataset

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
samples = pd.read_table(config["samples"], dtype = str, comment='#').set_index("sample", drop=False)
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

########## ANGSD ###########

# The following function is useful for now, but not robust to duplicate 
# names across column types. Needs improvement

def get_samples_from_pop(population):
    pop = population
    if pop == "all":
        return samples.index.values.tolist()
    elif pop == "all_excl_pca-admix":
        excl = config["excl_pca-admix"]
        return [s for s in samples.index.values.tolist() if s not in excl]
    elif pop in samples.depth.values:
        return samples.index[samples.depth == pop].values.tolist()
    elif pop in samples.population.values and pop not in samples.index:
        return samples.index[samples.population == pop].values.tolist()
    elif pop in samples.index and pop not in samples.population.values:
        return [pop]

def get_bamlist_bams(wildcards):
    pop = wildcards.population
    return expand(results+"/bams/{sample}{{dp}}.bam", 
                sample = get_samples_from_pop(pop))

def get_bamlist_bais(wildcards):
    pop = wildcards.population
    return expand(results+"/bams/{sample}{{dp}}.bam.bai", 
                sample = get_samples_from_pop(pop))