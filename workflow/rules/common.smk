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
    df = pd.DataFrame(contigs, columns=['contig','length']).set_index('contig')
    df.drop(index=config["reference"]["mito"]+
        config["reference"]["sex-linked"]+config["reference"]["exclude"],
        inplace = True)
    if chunk_size < max(df['length']):
        raise ValueError("Config invalid - chunk_size ("+str(chunk_size)+
            ") cannot be smaller than the largest contig in the reference ("+
            str(max(df['length']))+"). Please set chunk_size to a value "+
            "greater than or equal to "+str(max(df['length']))+").")
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
            if total > chunk_size:
                new_df = pd.DataFrame(included)
                dfs.append(new_df)
                included = []
                #new_df = df.iloc[0:0, :].copy()
                total = row['length']
            included.append(row)
            if n+1 == len(df):
                new_df = pd.DataFrame(included)
                dfs.append(new_df)    
    return dfs

chunks = chunkify(config["reference"]["fasta"], config["chunk_size"])
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

# Pre-processing

## Get fastq inputs for fastp
def get_raw_fastq(wildcards):
    unit = units.loc[wildcards.sample, ["fq1", "fq2"]]
    return {'r1': unit.fq1, 'r2': unit.fq2}

# Reference

## Get repeatmasker inputs
def get_repmaskin(wildcards):
	dic = {'ref' : "results/ref/{ref}/{ref}.fa"}
	if config["analyses"]["repeatmasker"]["local_lib"]:
		dic.update({'lib': config["analyses"]["repeatmasker"]["local_lib"]})
	elif config["analyses"]["repeatmasker"]["build_lib"]:
		dic.update({'lib': rules.repeatmodeler.output.fa})
	return dic

## Get inputs for combined filters file

def get_bed_filts(wildcards):
	bedin = []
	bedsum = []
	# add minimum size filter if set
	if config["reference"]["min_size"] and config["reference"]["min_size"] > 0:
		bedin.extend(expand("results/datasets/{{dataset}}/filters/small_scaffs/{{ref}}_scaff{size}bp.bed",size=config["reference"]["min_size"]))
		bedsum.extend(expand("results/datasets/{{dataset}}/filters/small_scaffs/{{ref}}_scaff{size}bp.bed.sum",size=config["reference"]["min_size"]))
	# add sex chromosome filter if set
	if config["reference"]["sex-linked"] or config["reference"]["exclude"] or config["reference"]["mito"]:
		bedin.append("results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_excl.bed")
		bedsum.append("results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_excl.bed.sum")
	# add mappability filter if set
	if config["analyses"]["genmap"]:
		bedin.append("results/datasets/{dataset}/filters/lowmap/{ref}_lowmap.bed")
		bedsum.append("results/datasets/{dataset}/filters/lowmap/{ref}_lowmap.bed.sum")
	# add repeat filter if set
	if any(config["analyses"]["repeatmasker"].values()):
		bedin.append("results/ref/{ref}/repeatmasker/{ref}.fa.out.gff")
		bedsum.append("results/ref/{ref}/repeatmasker/{ref}.fa.out.gff.sum")
	# add global depth extremes filter if set
	if config["analyses"]["extreme_depth"]:
		bedin.append("results/datasets/{dataset}/filters/depth/{dataset}.{ref}{dp}_extreme-depth.bed")
		bedsum.append("results/datasets/{dataset}/filters/depth/{dataset}.{ref}{dp}_extreme-depth.bed.sum")
	# add dataset level missing data filter if set
	if config["analyses"]["dataset_missing_data"]:
		bedin.extend(expand("results/datasets/{{dataset}}/filters/missdata/{{dataset}}.{{ref}}_all{{dp}}_under{miss}.bed",miss=config["analyses"]["dataset_missing_data"]))
		bedsum.extend(expand("results/datasets/{{dataset}}/filters/missdata/{{dataset}}.{{ref}}_all{{dp}}_under{miss}.bed.sum",miss=config["analyses"]["dataset_missing_data"]))
	# add population level missing data filter if set
	if config["analyses"]["population_missing_data"]:
		bedin.extend(expand("results/datasets/{{dataset}}/filters/missdata/{{dataset}}.{{ref}}_{population}{{dp}}_under{miss}.bed",
				miss=config["analyses"]["population_missing_data"],
				population=pop_list))
		bedsum.extend(expand("results/datasets/{{dataset}}/filters/missdata/{{dataset}}.{{ref}}_{population}{{dp}}_under{miss}.bed.sum",
				miss=config["analyses"]["population_missing_data"],
				population=pop_list))
	return {'gen': "results/ref/{ref}/beds/genome.bed", 
			'sum': "results/ref/{ref}/beds/genome.bed.sum", 
			'filt': bedin, 'sums': bedsum}

# Mapping

## Get read groups for mapping
def get_read_group(wildcards):
    return r"-R '@RG\tID:{unit}\tSM:{sample}\tLB:{sample}\tPL:{platform}'"\
        .format(
            unit=units.loc[wildcards.sample, "unit"],
            sample=wildcards.sample,
            platform=units.loc[wildcards.sample, "platform"]
        )

## Select which duplicate removal process the bam goes through
def get_dedup_bam(wildcards):
    s = wildcards.sample
    if s in samples.index[samples.time == "historical"]:
        return ["results/mapping/dedup/{sample}.{ref}.merged.rmdup.bam",
                "results/mapping/dedup/{sample}.{ref}.merged.rmdup.bam.bai"]
    elif s in samples.index[samples.time == "modern"]:
        return ["results/mapping/dedup/{sample}.{ref}.clipped.rmdup.bam",
                "results/mapping/dedup/{sample}.{ref}.clipped.rmdup.bam.bai"]

## Select final bam file to symlink
def get_final_bam(wildcards):
    # Determines if bam should use Picard or DeDup for duplicate removal
    s = wildcards.sample
    if s in samples.index[samples.time == "historical"]:
        return {"bam": "results/mapping/bams/{sample}.{ref}.rmdup.realn.rescaled.bam",
                "bai": "results/mapping/bams/{sample}.{ref}.rmdup.realn.rescaled.bam.bai"}
    elif s in samples.index[samples.time == "modern"]:
        return {"bam": "results/mapping/bams/{sample}.{ref}.rmdup.realn.bam",
                "bai": "results/mapping/bams/{sample}.{ref}.rmdup.realn.bam.bai"}

## Get flagstat files to estimate endogenous content
def get_endo_cont_stat(wildcards):
    # Gets bam file for endogenous content calculation
    s = wildcards.sample
    if s in samples.index[samples.time == "modern"]:
        return "results/mapping/mapped/{sample}.{ref}.paired.flagstat"
    elif s in samples.index[samples.time == "historical"]:
        return "results/mapping/mapped/{sample}.{ref}.merged.flagstat"

# ANGSD

## Get options for making beagle files, depends on whether the beagle file is 
## for the whole dataset (all filtered sites go in and SNPs are called) or for 
## a single population (SNPs previously found in whole dataset go in)
def get_snpset(wildcards):
	pop = wildcards.population
	if pop == "all":
		return ["results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.sites",
				"results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.sites.idx"]
	else:
		return ["results/datasets/{dataset}/filters/snps/{dataset}.{ref}{dp}_snps.sites",
				"results/datasets/{dataset}/filters/snps/{dataset}.{ref}{dp}_snps.sites.idx"]

def get_popopts(wildcards):
	pop = wildcards.population
	if pop == "all":
		return "-doMajorMinor 1 -SNP_pval "+ \
				str(config["params"]["angsd"]["snp_pval"])+" -minMaf "+ \
				str(config["params"]["angsd"]["min_maf"])
	else:
		return "-doMajorMinor 3"

## Remove requested individuals from beagle for pca and admix
def get_excl_ind_cols(wildcards):
	exclinds = config["excl_pca-admix"]
	exclindex = [samples.index.to_list().index(i) for i in exclinds]
	col1 = [x*3+4 for x in exclindex]
	col2 = [x+1 for x in col1]
	col3 = [x+1 for x in col2]
	remove = col1+col2+col3
	remove_string = ','.join([str(i) for i in remove])
	return remove_string

# Get all possible kinship estimate pairings
def get_kinship(wildcards):
	combos = list(itertools.combinations(samples.index, 2))
	# sort inds alphebetically, this ensures that should new inds be added
	# after generating some SFS, the reordering of the combinations won't
	# lead to generating identical SFS with the individuals swapped
	combos = [sorted(pair) for pair in combos]
	ind1 = [pair[0] for pair in combos]
	ind2 = [pair[1] for pair in combos]
	return expand("results/datasets/{{dataset}}/analyses/kinship/waples2019/{{dataset}}.{{ref}}_{ind1}-{ind2}{{dp}}.kinship",
				zip, ind1=ind1, ind2=ind2)

def get_total_bed(wildcards):
	if wildcards.prefix == "results/mapping/qc/ind_depth/unfiltered/":
		return "results/ref/{ref}/beds/genome.bed"
	elif wildcards.prefix == "results/datasets/{dataset}/qc/ind_depth/filtered/":
		return "results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.bed"

def get_depth_header(wildcards):
	if wildcards.prefix == "results/mapping/qc/ind_depth/unfiltered/":
		return "genome"
	elif wildcards.prefix == "results/datasets/{dataset}/qc/ind_depth/filtered/":
		return "filt"

def get_sample_qcs(wildcards):
	dic = {'inds': "results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
			'unfilt': "results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_all{dp}.depth.sum",
			'filt': "results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}.{ref}_all{dp}.depth.sum"}
	if config["analyses"]["endogenous_content"]:
		dic.update({'endo': "results/datasets/{dataset}/qc/endogenous_content/{dataset}.{ref}_all.endo.tsv"})
	return dic

def get_fst(wildcards):
	if wildcards.unit == "ind":
		unit = samples.index
	elif wildcards.unit == "pop":
		unit = pop_list
	combos = list(itertools.combinations(unit, 2))
	# sort pops alphebetically, this ensures that should new pops be added
	# after generating some SFS, the reordering of the combinations won't
	# lead to generating identical SFS with the populations swapped
	combos = [sorted(pair) for pair in combos]
	pop1 = [pair[0] for pair in combos]
	pop2 = [pair[1] for pair in combos]
	return expand("results/datasets/{{dataset}}/analyses/fst/{{dataset}}.{{ref}}_{population1}-{population2}{{dp}}.fst.global",
				zip, population1=pop1, population2=pop2)

def get_auto_sum(wildcards):
	if config["reference"]["sex-linked"] or config["reference"]["exclude"] or config["reference"]["mito"]:
		return "results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_excl.bed.sum"
	else:
		return "results/ref/{ref}/beds/genome.bed"

# list samples in a group
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

# list bam files for a grouping

def get_bamlist_bams(wildcards):
    pop = wildcards.population
    return expand("results/datasets/{{dataset}}/bams/{sample}.{{ref}}{{dp}}.bam", 
                sample = get_samples_from_pop(pop))

# list bai files for a grouping

def get_bamlist_bais(wildcards):
    pop = wildcards.population
    return expand("results/datasets/{{dataset}}/bams/{sample}.{{ref}}{{dp}}.bam.bai", 
                sample = get_samples_from_pop(pop))