# A set of common functions and variables that may be called in various Snakefiles


# Import required modules for workflow

from snakemake.utils import validate, min_version
import pandas as pd
import os
import itertools
import csv

# Check for minimum Snakemake version

min_version("7.25.0")

# Load default configfile


configfile: "config/config.yaml"


# Create variables for software containers (for easier version updating)

angsd_container = "docker://zjnolen/angsd:0.940"
pcangsd_container = "docker://zjnolen/pcangsd:1.10"
evaladmix_container = "docker://zjnolen/evaladmix:0.961"
ngsf_hmm_container = "docker://zjnolen/ngsf-hmm:1.1.0"
mapdamage_container = "docker://quay.io/biocontainers/mapdamage2:2.2.1--pyr40_0"
ngsrelate_container = "docker://zjnolen/ngsrelate:20220925-ec95c8f"
ngsld_container = "docker://zjnolen/ngsld:1.2.0"


# Define function for genome chunks to break up analysis (for parallelization)


def chunkify(reference_fasta, chunk_size):
    # Open reference fasta and get list of contigs and lengths
    contigs = []
    with open(reference_fasta, "r") as fasta_in:
        for header, seq in itertools.groupby(fasta_in, lambda x: x.startswith(">")):
            if header:
                contig = next(seq).strip(">").strip()
            seq_length = len("".join(seq).replace("\n", ""))
            if seq_length > 0:
                contigs = contigs + [[contig, seq_length]]
    df = pd.DataFrame(contigs, columns=["contig", "length"]).set_index("contig")
    # Drop excluded contigs from dataframe
    df.drop(
        index=config["reference"]["mito"]
        + config["reference"]["sex-linked"]
        + config["reference"]["exclude"],
        inplace=True,
    )
    # Set to 1 chunk if no chunk size provided
    if not chunk_size:
        chunk_size = sum(df["length"]) + 1
    # Error out if configured chunk size is smaller than the largest contig
    if chunk_size < max(df["length"]):
        raise ValueError(
            f"Config invalid - chunk_size ({str(chunk_size)}) cannot be smaller than "
            f"the largest contig in the reference ({max(df['length'])}). Please set "
            f"chunk_size to a value greater than or equal to {max(df['length'])}."
        )
    # Prepare list of chunk dataframes
    dfs = []
    included = []
    total = 0
    # Set minimum size to include contig
    if config["reference"]["min_size"]:
        minsize = config["reference"]["min_size"]
    else:
        minsize = 0
    # Add contigs to chunks, starting new chunk when adding a contig would take it over
    # the configured chunk size
    for n, (index, row) in enumerate(df.iterrows()):
        if row["length"] >= minsize:
            total += row["length"]
            if total > chunk_size:
                new_df = pd.DataFrame(included)
                dfs.append(new_df)
                included = []
                # new_df = df.iloc[0:0, :].copy()
                total = row["length"]
            included.append(row)
        if ((n + 1) == len(df)) and included:
            new_df = pd.DataFrame(included)
            dfs.append(new_df)
    return dfs


# Set up chunk variables

chunks = chunkify(config["reference"]["fasta"], config["chunk_size"])
chunklist = list(range(1, len(chunks) + 1))

# Load sample sheet

samples = pd.read_table(config["samples"], dtype=str, comment="#").set_index(
    "sample", drop=False
)


# Drop samples specified in config

samples = samples.drop(config["exclude_ind"])


# Load unit sheet and drop any unused units

units = pd.read_table(config["units"])
units = units[units["sample"].isin(list(samples.index))]


# Get lists of samples where bams are provided by users and where they are made
# by the pipeline

if "bam" in units:
    userbams = list(units["sample"][units["bam"].notnull()].unique())
    pipebams = list(units["sample"][units["bam"].isnull()].unique())


# Get a list of all the populations

if config["populations"]:
    pop_list = config["populations"]
else:
    pop_list = samples.population.unique().tolist()


# Get a list of user defined filter sets if present

filters = []
if not config["only_filter_beds"]:
    filters = filters + ["allsites"]
if any(config["filter_beds"].values()):
    filters = filters + list(config["filter_beds"].keys())


# Define various helper functions for Snakefiles


# Pre-processing


## Get fastq inputs for fastp
def get_raw_fastq(wildcards):
    unit = units.loc[
        (units["sample"] == wildcards.sample)
        & (units["unit"] == wildcards.unit)
        & (units["lib"] == wildcards.lib),
        ["sample", "fq1", "fq2"],
    ].set_index("sample")
    return {"sample": [unit.fq1[0], unit.fq2[0]]}


## Get minimum overlap to collapse read pairs per sample
def get_min_overlap(wildcards):
    s = wildcards.sample
    if s in samples.index[samples.time == "historical"]:
        return config["params"]["fastp"]["min_overlap_hist"]
    elif s in samples.index[samples.time == "modern"]:
        return config["params"]["fastp"]["min_overlap_mod"]


# Reference


## Get repeatmasker inputs
def get_repmaskin(wildcards):
    dic = {"ref": "results/ref/{ref}/{ref}.fa"}
    if config["analyses"]["repeatmasker"]["local_lib"]:
        dic.update({"lib": config["analyses"]["repeatmasker"]["local_lib"]})
    elif config["analyses"]["repeatmasker"]["build_lib"]:
        dic.update({"lib": rules.repeatmodeler.output.fa})
    return dic


## Get inputs for combined filters file
def get_bed_filts(wildcards):
    bedin = []
    bedsum = []
    # add minimum size filter if set
    if config["reference"]["min_size"] and config["reference"]["min_size"] > 0:
        bedin.extend(
            expand(
                "results/datasets/{{dataset}}/filters/small_scaffs/{{ref}}_scaff{size}bp.bed",
                size=config["reference"]["min_size"],
            )
        )
        bedsum.extend(
            expand(
                "results/datasets/{{dataset}}/filters/small_scaffs/{{ref}}_scaff{size}bp.bed.sum",
                size=config["reference"]["min_size"],
            )
        )
    # add sex chromosome filter if set
    if (
        config["reference"]["sex-linked"]
        or config["reference"]["exclude"]
        or config["reference"]["mito"]
    ):
        bedin.append(
            "results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_excl.bed"
        )
        bedsum.append(
            "results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_excl.bed.sum"
        )
    # add mappability filter if set
    if config["analyses"]["pileup-mappability"]:
        bedin.extend(
            expand(
                "results/datasets/{{dataset}}/filters/pileupmap/{{ref}}_k{k}_e{e}_{thresh}.bed",
                k=config["params"]["genmap"]["K"],
                e=config["params"]["genmap"]["E"],
                thresh=config["params"]["genmap"]["map_thresh"],
            )
        )
        bedsum.extend(
            expand(
                "results/datasets/{{dataset}}/filters/pileupmap/{{ref}}_k{k}_e{e}_{thresh}.bed.sum",
                k=config["params"]["genmap"]["K"],
                e=config["params"]["genmap"]["E"],
                thresh=config["params"]["genmap"]["map_thresh"],
            )
        )
    # add repeat filter if set
    if any(config["analyses"]["repeatmasker"].values()):
        bedin.append("results/ref/{ref}/repeatmasker/{ref}.fa.out.gff")
        bedsum.append("results/ref/{ref}/repeatmasker/{ref}.fa.out.gff.sum")
    # add global depth extremes filter if set
    if config["analyses"]["extreme_depth"]:
        if (
            not config["params"]["extreme_depth_filt"]["filt-on-dataset"]
            and not config["params"]["extreme_depth_filt"]["filt-on-depth-classes"]
        ):
            raise ValueError(
                f"Config invalid - extreme_depth filter is set to true, but neither "
                f"filt-on-dataset or filt-on-depth-classes is set to true in the "
                f"params, so no groups are defined for the filter to act on. Please "
                f"set the filter to false if you do not want to filter on depth. If "
                f"you do, please set the categories the filter will happen on."
            )
        if config["params"]["extreme_depth_filt"]["filt-on-dataset"]:
            bedin.append(
                "results/datasets/{dataset}/filters/depth/{dataset}.{ref}_all{dp}_extreme-depth.bed"
            )
            bedsum.append(
                "results/datasets/{dataset}/filters/depth/{dataset}.{ref}_all{dp}_extreme-depth.bed.sum"
            )
        if config["params"]["extreme_depth_filt"]["filt-on-depth-classes"]:
            bedin.extend(
                expand(
                    "results/datasets/{{dataset}}/filters/depth/{{dataset}}.{{ref}}_{population}{{dp}}_extreme-depth.bed",
                    population=list(set(samples.depth.values)),
                )
            )
            bedsum.extend(
                expand(
                    "results/datasets/{{dataset}}/filters/depth/{{dataset}}.{{ref}}_{population}{{dp}}_extreme-depth.bed.sum",
                    population=list(set(samples.depth.values)),
                )
            )
    # add dataset level missing data filter if set
    if config["analyses"]["dataset_missing_data"]:
        bedin.extend(
            expand(
                "results/datasets/{{dataset}}/filters/missdata/{{dataset}}.{{ref}}_all{{dp}}_under{miss}.bed",
                miss=config["analyses"]["dataset_missing_data"],
            )
        )
        bedsum.extend(
            expand(
                "results/datasets/{{dataset}}/filters/missdata/{{dataset}}.{{ref}}_all{{dp}}_under{miss}.bed.sum",
                miss=config["analyses"]["dataset_missing_data"],
            )
        )
    # add population level missing data filter if set
    if config["analyses"]["population_missing_data"]:
        bedin.extend(
            expand(
                "results/datasets/{{dataset}}/filters/missdata/{{dataset}}.{{ref}}_{population}{{dp}}_under{miss}.bed",
                miss=config["analyses"]["population_missing_data"],
                population=pop_list,
            )
        )
        bedsum.extend(
            expand(
                "results/datasets/{{dataset}}/filters/missdata/{{dataset}}.{{ref}}_{population}{{dp}}_under{miss}.bed.sum",
                miss=config["analyses"]["population_missing_data"],
                population=pop_list,
            )
        )
    return {
        "gen": "results/ref/{ref}/beds/genome.bed",
        "sum": "results/ref/{ref}/beds/genome.bed.sum",
        "filt": bedin,
        "sums": bedsum,
    }


# Get bed file for user defined filters
def get_newfilt(wildcards):
    return config["filter_beds"][wildcards.sites]


# Mapping


## Get read groups for mapping
def get_read_group(wildcards):
    return r"-R '@RG\tID:{id}\tPU:{pu}\tSM:{sample}\tLB:{lib}\tPL:{platform}'".format(
        id=f"{wildcards.unit}.{wildcards.lib}",
        pu=f"{wildcards.unit}",
        sample=wildcards.sample,
        lib=wildcards.lib,
        platform=units.loc[
            (units["sample"] == wildcards.sample)
            & (units["unit"] == wildcards.unit)
            & (units["lib"] == wildcards.lib),
            "platform",
        ].to_list(),
    )


## Get single unit/lib bams for merging
def get_sample_bams(wildcards):
    reads = units.loc[units["sample"] == wildcards.sample]
    combos = reads[["sample", "unit", "lib"]].agg("_".join, axis=1)
    return expand(
        "results/mapping/mapped/{combo}.{{ref}}.{{pairing}}.bam", combo=combos
    )


## Select which duplicate removal process the bam goes through
def get_dedup_bams(wildcards):
    if config["analyses"]["mapping"]["historical_only_collapsed"]:
        s = wildcards.sample
        if s in samples.index[samples.time == "historical"]:
            return ["results/mapping/dedup/{sample}.{ref}.merged.rmdup.bam"]
        elif s in samples.index[samples.time == "modern"]:
            return [
                "results/mapping/dedup/{sample}.{ref}.paired.rmdup.bam",
                "results/mapping/dedup/{sample}.{ref}.merged.rmdup.bam",
            ]
    else:
        return [
            "results/mapping/dedup/{sample}.{ref}.paired.rmdup.bam",
            "results/mapping/dedup/{sample}.{ref}.merged.rmdup.bam",
        ]


## Function to pull user provided bam. For now, error raising is handled in
## get_final_bam, as that runs before this.
def get_user_bam(wildcards):
    s = wildcards.sample
    return {"bam": units.loc[units["sample"] == s, "bam"].values[0]}


## Determine what bam to use in analyses. This decides whether to use user provided
## bams, and process raw reads if not. It also determines if mapdamage rescaling will
## be performed on historical bams.
def get_final_bam(wildcards):
    s = wildcards.sample
    if ("bam" in units) and any(units.loc[units["sample"] == s, "bam"].notnull()):
        if sum(units["sample"] == s) > 1:
            raise ValueError(
                f"Config invalid - sample {s} appears in units file multiple times, "
                f"but at least one of these times has a user provided bam. Multiple "
                f"units are only supported when the pipeline is doing bam processing "
                f"for a sample. Please process the sample into one bam file, or "
                f"process all the raw reads using the pipeline without supplying a bam "
                f"for the sample."
            )
        if config["params"]["clipoverlap"]["clip_user_provided_bams"]:
            return {
                "bam": "results/datasets/{dataset}/bams/clipped_user_bams/{sample}.{ref}.clip.bam",
                "bai": "results/datasets/{dataset}/bams/clipped_user_bams/{sample}.{ref}.clip.bam.bai",
            }
        return {
            "bam": "results/datasets/{dataset}/bams/user_bams/{sample}.{ref}.user-processed.bam",
            "bai": "results/datasets/{dataset}/bams/user_bams/{sample}.{ref}.user-processed.bam.bai",
        }
    if (s in samples.index[samples.time == "historical"]) and (
        config["analyses"]["mapdamage_rescale"]
    ):
        return {
            "bam": "results/mapping/bams/{sample}.{ref}.rmdup.realn.rescaled.clip.bam",
            "bai": "results/mapping/bams/{sample}.{ref}.rmdup.realn.rescaled.clip.bam.bai",
        }
    return {
        "bam": "results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.bam",
        "bai": "results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.bam.bai",
    }


## Get flagstat file for endogenous content calculation
def get_endo_cont_stat(wildcards):
    s = wildcards.sample
    if userbams and s in userbams:
        return {
            "paired": "results/datasets/{dataset}/bams/user_bams/{sample}.{ref}.user-processed.flagstat",
            "merged": "results/datasets/{dataset}/bams/user_bams/{sample}.{ref}.user-processed.flagstat"
        }
    if (config["analyses"]["mapping"]["historical_only_collapsed"]) and (
        s in samples.index[samples.time == "historical"]
    ):
        return {
            "paired": "results/mapping/mapped/{sample}.{ref}.merged.flagstat",
            "merged": "results/mapping/mapped/{sample}.{ref}.merged.flagstat",
        }
    else:
        return {
            "paired": "results/mapping/mapped/{sample}.{ref}.paired.flagstat",
            "merged": "results/mapping/mapped/{sample}.{ref}.merged.flagstat",
        }


# ANGSD

## Get GLF input for beagle and SAF. If statistics are only desired for sites from a
## user defined BED file, GLFs will only be made for those sites to prevent estimating
## GLs for entire genome when only a subset is desired.


# def get_glf(wildcards):
#     if config["only_filter_beds"] and not any(config["filter_beds"].values()):
#         raise ValueError(
#             f"Config invalid - 'only_filter_beds' cannot be true without supplying bed "
#             f"files to at least one 'filter_beds' key."
#         )
#     elif config["only_filter_beds"]:
#         return "results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.glf.gz"
#     else:
#         return "results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_allsites-filts.glf.gz"


## Get random sampling proportion depending on if LD decay is being calculated
## or if LD pruning is being done
def get_ngsld_sampling(wildcards):
    if wildcards.path == "beagles/pruned/ngsLD":
        return "1"
    elif wildcards.path == "analyses/ngsLD/chunks":
        return config["params"]["ngsld"]["rnd_sample"]


## Get sample size for r^2 sample size corrections on LD decay
def get_ngsld_n(wildcards):
    if config["params"]["ngsld"]["fit_LDdecay_n_correction"]:
        nind = len(get_samples_from_pop(wildcards.population))
        return f"--n_ind {nind}"
    else:
        return ""


## Remove requested individuals from beagle for pca and admix
def get_excl_ind_cols(wildcards):
    exclinds = config["excl_pca-admix"]
    exclindex = [samples.index.to_list().index(i) for i in exclinds]
    col1 = [x * 3 + 4 for x in exclindex]
    col2 = [x + 1 for x in col1]
    col3 = [x + 1 for x in col2]
    remove = col1 + col2 + col3
    remove_string = ",".join([str(i) for i in remove])
    return remove_string


## Get all possible kinship estimate pairings
def get_kinship(wildcards):
    combos = list(itertools.combinations(samples.index, 2))
    # sort inds alphebetically, this ensures that should new inds be added
    # after generating some SFS, the reordering of the combinations won't
    # lead to generating identical SFS with the individuals swapped
    combos = [sorted(pair) for pair in combos]
    ind1 = [pair[0] for pair in combos]
    ind2 = [pair[1] for pair in combos]
    return expand(
        "results/datasets/{{dataset}}/analyses/kinship/ibsrelate_sfs/{{dataset}}.{{ref}}_{ind1}-{ind2}{{dp}}_{{sites}}-filts.kinship",
        zip,
        ind1=ind1,
        ind2=ind2,
    )


## Get bedfile for whole genome or filtered genome, to set the denomintor of coverage
## calculation
def get_total_bed(wildcards):
    if wildcards.prefix == "results/mapping/qc/ind_depth/unfiltered/":
        return "results/ref/{ref}/beds/genome.bed"
    elif (
        wildcards.prefix
        == f"results/datasets/{wildcards.dataset}/qc/ind_depth/filtered/"
    ):
        return "results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{group}.bed"


## Gather sample QC output files for concatenation
def get_sample_qcs(wildcards):
    dic = {
        "inds": "results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
        "unfilt": "results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_all{dp}_allsites-unfilt.depth.sum",
        "filt": expand(
            "results/datasets/{{dataset}}/qc/ind_depth/filtered/{{dataset}}.{{ref}}_all{{dp}}_{sites}-filts.depth.sum",
            sites=filters,
        ),
        "endo": "results/datasets/{dataset}/qc/endogenous_content/{dataset}.{ref}_all.endo.tsv",
    }
    return dic


## Get a list of all pairwise combos of a set of items
def pairwise_combos(items):
    combos = list(itertools.combinations(items, 2))
    # sort pops alphebetically, this ensures that should new pops be added
    # after generating some SFS, the reordering of the combinations won't
    # lead to generating identical SFS with the populations swapped
    combos = [sorted(pair) for pair in combos]
    return combos


## Get fst tables for population aggregation, returning all pairwise combinations of
## individuals or populations, for global or windowed estimates
def get_fst(wildcards):
    if wildcards.unit == "ind":
        unit = samples.index
    elif wildcards.unit == "pop":
        unit = pop_list
    combos = pairwise_combos(unit)
    pop1 = [pair[0] for pair in combos]
    pop2 = [pair[1] for pair in combos]
    if wildcards.scale == "global":
        return expand(
            "results/datasets/{{dataset}}/analyses/fst/{{dataset}}.{{ref}}_{population1}-{population2}{{dp}}_{{sites}}-filts.fst.global.tsv",
            zip,
            population1=pop1,
            population2=pop2,
        )
    elif wildcards.scale == "window":
        return expand(
            "results/datasets/{{dataset}}/analyses/fst/{{dataset}}.{{ref}}_{population1}-{population2}{{dp}}_{{sites}}-filts.fst.window_{{win}}_{{step}}.tsv",
            zip,
            population1=pop1,
            population2=pop2,
        )


## Get summary file that describes the autosomes (depends on if non-autosomal contigs
## have been specified in the config)
def get_auto_sum(wildcards):
    if (
        config["reference"]["sex-linked"]
        or config["reference"]["exclude"]
        or config["reference"]["mito"]
    ):
        return (
            "results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_excl.bed.sum"
        )
    else:
        return "results/ref/{ref}/beds/genome.bed"


# List samples in a group
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


# List bam files for a grouping
def get_bamlist_bams(wildcards):
    pop = wildcards.population
    return expand(
        "results/datasets/{{dataset}}/bams/{sample}.{{ref}}{{dp}}.bam",
        sample=get_samples_from_pop(pop),
    )


# List bai files for a grouping
def get_bamlist_bais(wildcards):
    pop = wildcards.population
    return expand(
        "results/datasets/{{dataset}}/bams/{sample}.{{ref}}{{dp}}.bam.bai",
        sample=get_samples_from_pop(pop),
    )


# Get string to describe depth subsampling for report
def dp_report(wildcards):
    dp = wildcards.dp
    if dp == "":
        return {"Subsampling": "None"}
    else:
        return {"Subsampling": f"{dp.replace('.dp','')}X"}


# Get string to describe units for Fst for report
def unit_report(wildcards):
    unit = wildcards.unit
    if unit == "ind":
        return {"Unit": "Individuals"}
    elif unit == "pop":
        return {"Unit": "Populations"}


# Get string to describe summary statistic for thetas for report
def theta_report(wildcards):
    stat = wildcards.stat
    if stat == "watterson":
        return "Watterson's Theta"
    elif stat == "pi":
        return "Nucleotide Diversity (Pi)"
    elif stat == "tajima":
        return "Tajima's D"
