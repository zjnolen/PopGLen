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

angsd_container = "docker://quay.io/biocontainers/angsd:0.940--hf5e1c6e_3"
bamutil_container = "docker://quay.io/biocontainers/bamutil:1.0.15--h43eeafb_5"
bedops_container = "docker://quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0"
bedtools_container = "docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2"
damageprofiler_container = (
    "docker://quay.io/biocontainers/damageprofiler:1.1--hdfd78af_2"
)
dedup_container = "docker://quay.io/biocontainers/dedup:0.12.8--hdfd78af_1"
evaladmix_container = "docker://ghcr.io/zjnolen/evaladmix:0.961"
genmap_container = "docker://quay.io/biocontainers/genmap:1.3.0--he1b5a44_0"
multiqc_container = "docker://quay.io/biocontainers/multiqc:1.23--pyhdfd78af_0"
ngsf_hmm_container = "docker://ghcr.io/zjnolen/ngsf-hmm:1.1.0"
ngsld_container = "docker://ghcr.io/zjnolen/ngsld:1.2.0"
ngsrelate_container = "docker://ghcr.io/zjnolen/ngsrelate:20220925-ec95c8f"
pandas_container = "docker://quay.io/biocontainers/pandas:2.2.1"
pcangsd_container = "docker://ghcr.io/zjnolen/pcangsd:1.35"
qualimap_container = "docker://quay.io/biocontainers/qualimap:2.3--hdfd78af_0"
r_container = "docker://ghcr.io/zjnolen/popglen-r:develop-241016"
repeatmodmask_container = (
    "docker://quay.io/biocontainers/repeatmodeler:2.0.5--pl5321hdfd78af_0"
)
samtools_container = "docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_1"
shell_container = "docker://ghcr.io/zjnolen/popglen-shell:0.4.0"

# Define function for genome chunks to break up analysis (for parallelization)


def chunkify(reference_fasta, chunk_size):
    # Open reference fasta and get list of contigs and lengths
    contigs = []
    with open(reference_fasta, "r") as fasta_in:
        for header, seq in itertools.groupby(fasta_in, lambda x: x.startswith(">")):
            if header:
                contig = next(seq).strip(">").strip().split(maxsplit=1)[0]
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
else:
    userbams = []
    pipebams = list(samples.index)


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


###############################################################################
# ========== Define various helper functions for Snakefiles ================= #
###############################################################################

# Pre-processing


## Get fastq inputs for fastp
def get_raw_fastq(wildcards):
    if "fq1" in units:
        unit = units.loc[
            (units["sample"] == wildcards.sample)
            & (units["unit"] == wildcards.unit)
            & (units["lib"] == wildcards.lib),
            ["sample", "fq1", "fq2"],
        ].set_index("sample")
        if not pd.isna(unit.fq1.item()):
            return {"sample": [unit.fq1.item(), unit.fq2.item()]}
    if "sra" in units:
        sra = units[
            (units["sample"] == wildcards.sample)
            & (units["unit"] == wildcards.unit)
            & (units["lib"] == wildcards.lib)
        ].sra.item()
        if not pd.isna(sra):
            return {
                "sample": [
                    f"results/downloaded_fastq/{sra}_1.fastq.gz",
                    f"results/downloaded_fastq/{sra}_2.fastq.gz",
                ]
            }


## Get correct reports to compile a fastp multiqc report
def multiqc_input_fastp(wildcards):
    reports = []
    # Check if pipeline is actually processing any fastq files
    if len(pipebams) > 0:
        # subset units to samples that are starting at fastq
        pipeunits = units[units["sample"].isin(pipebams)]
        # join with sample list to know 'historical' or 'modern' sample context
        pipeunits = pd.merge(pipeunits, samples, left_on="sample", right_index=True)
        # add historical, merged fastq to report
        histunits = pipeunits[pipeunits["time"] == "historical"]
        reports.extend(
            expand(
                "results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_merged.json",
                zip,
                sample=histunits["sample"].tolist(),
                unit=histunits["unit"].tolist(),
                lib=histunits["lib"].tolist(),
            )
        )
        # add modern, paired fastq to report
        modunits = pipeunits[pipeunits["time"] == "modern"]
        reports.extend(
            expand(
                "results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_paired.json",
                zip,
                sample=modunits["sample"].tolist(),
                unit=modunits["unit"].tolist(),
                lib=modunits["lib"].tolist(),
            )
        )
    return reports


# Reference


## Get repeatmasker inputs
def get_repmaskin(wildcards):
    dic = {"ref": "results/ref/{ref}/{ref}.fa"}
    if config["analyses"]["repeatmasker"]["local_lib"]:
        dic.update({"lib": config["analyses"]["repeatmasker"]["local_lib"]})
    elif config["analyses"]["repeatmasker"]["build_lib"]:
        dic.update({"lib": rules.repeatmodeler.output.fa})
    return dic


## Determine if making repeat bed or using user supplied
def get_rep_file(wildcards):
    if config["analyses"]["repeatmasker"]["bed"]:
        return {"rep": config["analyses"]["repeatmasker"]["bed"]}
    return {"rep": "results/ref/{ref}/repeatmasker/{ref}.fa.out.gff"}


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
        bedin.append("results/ref/{ref}/repeatmasker/{ref}.fa.out.bed")
        bedsum.append("results/ref/{ref}/repeatmasker/{ref}.fa.out.sum")
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
    return r"'@RG\tID:{id}\tPU:{pu}\tSM:{sample}\tLB:{lib}\tPL:{platform}'".format(
        id=f"{wildcards.unit}.{wildcards.lib}",
        pu=f"{wildcards.unit}",
        sample=wildcards.sample,
        lib=wildcards.lib,
        platform=units.loc[
            (units["sample"] == wildcards.sample)
            & (units["unit"] == wildcards.unit)
            & (units["lib"] == wildcards.lib),
            "platform",
        ].to_list()[0],
    )


## Get single unit/lib bams for merging
def get_lib_bams(wildcards):
    reads = units.loc[units["sample"] == wildcards.sample]
    reads = reads.loc[reads["lib"] == wildcards.lib]
    combos = reads[["sample", "unit", "lib"]].agg("_".join, axis=1)
    return expand(
        "results/mapping/mapped/{combo}.{{ref}}.{aligner}.merged.bam",
        combo=combos,
        aligner=config["analyses"]["mapping"]["historical_collapsed_aligner"],
    )


def get_paired_bams(wildcards):
    reads = units.loc[units["sample"] == wildcards.sample]
    combos = reads[["sample", "unit", "lib"]].agg("_".join, axis=1)
    if wildcards.sample in samples.index[samples.time == "historical"]:
        return expand(
            "results/mapping/mapped/{combo}.{{ref}}.mem.uncollapsed.bam",
            combo=combos,
        )
    return expand(
        "results/mapping/mapped/{combo}.{{ref}}.mem.paired.bam",
        combo=combos,
    )


## Select which duplicate removal process the bam goes through
def get_dedup_bams(wildcards):
    s = wildcards.sample
    if s in samples.index[samples.time == "historical"]:
        libs = list(set(units.loc[units["sample"] == wildcards.sample]["lib"]))
        if config["analyses"]["mapping"]["historical_only_collapsed"]:
            return expand(
                "results/mapping/dedup/{{sample}}_{lib}.{{ref}}.merged.rmdup.bam",
                lib=libs,
            )
        output = ["results/mapping/dedup/{sample}.{ref}.paired.rmdup.bam"]
        output.extend(
            expand(
                "results/mapping/dedup/{{sample}}_{lib}.{{ref}}.merged.rmdup.bam",
                lib=libs,
            )
        )
        return output
    return "results/mapping/dedup/{sample}.{ref}.paired.rmdup.bam"


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
                "bam": "results/mapping/user-provided-bams/{sample}.{ref}.clip.bam",
                "bai": "results/mapping/user-provided-bams/{sample}.{ref}.clip.bam.bai",
            }
        return {
            "bam": "results/mapping/user-provided-bams/{sample}.{ref}.user-processed.bam",
            "bai": "results/mapping/user-provided-bams/{sample}.{ref}.user-processed.bam.bai",
        }
    if (s in samples.index[samples.time == "historical"]) and (
        config["analyses"]["mapdamage_rescale"]
    ):
        return {
            "bam": "results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.rescaled.bam",
            "bai": "results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.rescaled.bam.bai",
        }
    return {
        "bam": "results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.bam",
        "bai": "results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.bam.bai",
    }


# Sample QC


## Get flagstat file for endogenous content calculation
def get_endo_cont_stat(wildcards):
    s = wildcards.sample
    if userbams and s in userbams:
        return {
            "paired": "results/mapping/user-provided-bams/{sample}.{ref}.user-processed.flagstat",
            "merged": "results/mapping/user-provided-bams/{sample}.{ref}.user-processed.flagstat",
        }
    if s in samples.index[samples.time == "historical"]:
        libs = list(set(units.loc[units["sample"] == wildcards.sample]["lib"]))
        if config["analyses"]["mapping"]["historical_only_collapsed"]:
            return {
                "merged": expand(
                    "results/mapping/mapped/{{sample}}_{lib}.{{ref}}.merged.flagstat",
                    lib=libs,
                )
            }
        return {
            "paired": "results/mapping/mapped/{sample}.{ref}.paired.flagstat",
            "merged": expand(
                "results/mapping/mapped/{{sample}}_{lib}.{{ref}}.merged.flagstat",
                lib=libs,
            ),
        }
    return {"paired": "results/mapping/mapped/{sample}.{ref}.paired.flagstat"}


## Get bedfile for whole genome or filtered genome, to set the denomintor of coverage
## calculation
def get_filter_sum(wildcards):
    if (
        wildcards.prefix
        == f"results/datasets/{wildcards.dataset}/qc/ind_depth/filtered/"
    ):
        if (config["subsample_by"] != "sitefilt") and config["redo_depth_filts"]:
            return "results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{group}.sum"
        return "results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{group}.sum"
    return "results/ref/{ref}/beds/genome.bed.sum"


## Gather sample QC output files for concatenation
def get_sample_qcs(wildcards):
    dic = {
        "inds": "results/datasets/{dataset}/poplists/{dataset}_all{dp}.indiv.list",
        "endo": "results/datasets/{dataset}/qc/endogenous_content/{dataset}.{ref}_all{dp}.endo.tsv",
        "unfilt": "results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_all{dp}_allsites-unfilt.depth.sum",
        "mapqbaseqfilt": expand(
            "results/mapping/qc/ind_depth/mapq-baseq-filtered/{{dataset}}.{{ref}}_all{{dp}}_allsites-mapq{mapq}-baseq{baseq}-filt.depth.sum",
            mapq=config["mapQ"],
            baseq=config["baseQ"],
        ),
        "filt": expand(
            "results/datasets/{{dataset}}/qc/ind_depth/filtered/{{dataset}}.{{ref}}_all{{dp}}_{sites}-filts.depth.sum",
            sites=filters,
        ),
    }
    return dic


def multiqc_input_qualimap(wildcards):
    reports = []
    if len(pipebams) > 0:
        reports.extend(
            expand(
                "results/mapping/qc/qualimap/{sample}.{{ref}}",
                sample=pipebams,
            )
        )
    if len(userbams) > 0:
        reports.extend(
            expand(
                "results/datasets/{{dataset}}/qc/user-provided-bams/qualimap/{sample}.{{ref}}",
                sample=userbams,
            )
        )
    return reports


# DNA Damage


def multiqc_input_dnadmg(wildcards):
    reports = []
    if config["analyses"]["damageprofiler"]:
        if len(pipebams) > 0:
            if config["params"]["damageprofiler"]["profile_modern"]:
                samset = set(samples.index)
            else:
                samset = set(samples.index[samples["time"] == "historical"])
            reports.extend(
                expand(
                    "results/mapping/qc/damageprofiler/{histsample}.{{ref}}/dmgprof.json",
                    histsample=list(set(pipebams) & samset),
                )
            )
        if len(userbams) > 0:
            if config["params"]["damageprofiler"]["profile_modern"]:
                samset = set(samples.index)
            else:
                samset = set(samples.index[samples["time"] == "historical"])
            reports.extend(
                expand(
                    "results/mapping/qc/damageprofiler/{histsample}.{{ref}}.user-processed/dmgprof.json",
                    histsample=list(set(userbams) & samset),
                )
            )
    if config["analyses"]["mapdamage_rescale"]:
        if len(pipebams) > 0:
            reports.extend(
                expand(
                    "results/mapping/qc/mapdamage/{histsample}.{{ref}}/lgdistribution.txt",
                    histsample=list(
                        set(pipebams)
                        & set(samples.index[samples["time"] == "historical"])
                    ),
                )
            )
    return reports


# ANGSD


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


# Return sample list for making popfiles, ensuring dropped samples are left out
# in subsampled runs
def get_popfile_inds(wildcards):
    sam = get_samples_from_pop(wildcards.population)
    if wildcards.dp != "":
        sam = [s for s in sam if s not in config["drop_samples"]]
    return sam


# Return individual counts, dropping samples when required for subsampled runs
def get_nind(wildcards):
    sam = get_popfile_inds(wildcards)
    return len(sam)


# List bam files for a grouping
def get_bamlist_bams(wildcards):
    sam = get_popfile_inds(wildcards)
    return expand(
        "results/datasets/{{dataset}}/bams/{sample}.{{ref}}{{dp}}.bam",
        sample=sam,
    )


# List bai files for a grouping
def get_bamlist_bais(wildcards):
    sam = get_popfile_inds(wildcards)
    return expand(
        "results/datasets/{{dataset}}/bams/{sample}.{{ref}}{{dp}}.bam.bai",
        sample=sam,
    )


# Select filter file based off of full depth samples or subsampled depth
def filt_depth(wildcards):
    if (config["subsample_by"] != "sitefilt") and config["redo_depth_filts"]:
        return {
            "sites": "results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.sites",
            "idx": "results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.sites.idx",
        }
    return {
        "sites": "results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites",
        "idx": "results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites.idx",
    }


# Select mean depth to use for subsampling based on user preference. Sitefilt
# is recommended when trying to make sequencing depth uniform across samples in
# downstream analyses.
def depth_file(wildcards):
    if config["subsample_by"] == "unfilt":
        return "results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_{sample}_allsites-unfilt.depth.sum"
    if config["subsample_by"] == "mapqbaseq":
        return expand(
            "results/mapping/qc/ind_depth/mapq-baseq-filtered/{{dataset}}.{{ref}}_{{sample}}_allsites-mapq{mapq}-baseq{baseq}-filt.depth.sum",
            mapq=config["mapQ"],
            baseq=config["baseQ"],
        )
    if config["subsample_by"] == "sitefilt":
        return "results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}.{ref}_{sample}_allsites-filts.depth.sum"


# Choose to use an ancestral reference if present, otherwise use main reference
def get_anc_ref(wildcards):
    if config["ancestral"]:
        return {
            "anc": "results/ref/{ref}/{ref}.anc.fa",
            "ancfai": "results/ref/{ref}/{ref}.anc.fa.fai",
        }
    return {
        "anc": "results/ref/{ref}/{ref}.fa",
        "ancfai": "results/ref/{ref}/{ref}.fa.fai",
    }


# Define if minimum individual filtering will happen in ANGSD, and set the
# command line option
def get_minind(wildcards):
    pop = wildcards.population
    sam = get_popfile_inds(wildcards)
    if pop == "all":
        minind = int(len(sam) * config["params"]["angsd"]["minind_dataset"])
    else:
        minind = int(len(sam) * config["params"]["angsd"]["minind_pop"])
    minindopt = f"-minInd {minind}"
    return minindopt


# Determine if docounts is needed for beagle/maf calculation to keep it from
# slowing things down when it is not. It is only needed if the major and minor
# alleles are being inferred from counts (-doMajorMinor 2). This would also be
# needed if using -doMaf 8, but that is not supported currently by the pipeline.
def get_docounts(wildcards):
    if str(config["params"]["angsd"]["domajorminor"]) == "2":
        return "-doCounts 1"
    return ""


# Determine whether transitions should be removed based on user configuration
def get_trans(wildcards):
    if config["params"]["angsd"]["rmtrans"]:
        return 1
    return 0


# ngsLD


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
        sam = get_popfile_inds(wildcards)
        nind = len(sam)
        return "--n_ind %s" % nind
    else:
        return ""


# PCA/Admix


## Remove requested individuals from beagle for pca and admix
def get_excl_ind_cols(wildcards):
    exclinds = config["excl_pca-admix"]
    sam = get_popfile_inds(wildcards)
    exclinds = [s for s in exclinds if s in sam]
    exclindex = [sam.index(i) for i in exclinds]
    col1 = [x * 3 + 4 for x in exclindex]
    col2 = [x + 1 for x in col1]
    col3 = [x + 1 for x in col2]
    remove = col1 + col2 + col3
    remove_string = ",".join([str(i) for i in remove])
    return remove_string


# Kinship


## Get all possible kinship estimate pairings
def get_kinship(wildcards):
    sam = get_samples_from_pop("all")
    if wildcards.dp != "":
        sam = [s for s in sam if s not in config["drop_samples"]]
    combos = list(itertools.combinations(sam, 2))
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


# Fst


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
        sam = samples.index.tolist()
        if wildcards.dp != "":
            sam = [s for s in sam if s not in config["drop_samples"]]
        unit = sam
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


# Reports


# Get string to describe depth subsampling for report
def dp_report(wildcards):
    dp = wildcards.dp
    if dp == "":
        return {"Subsampling": "None"}
    else:
        return {"Subsampling": f"{dp.replace('.dp', '')}X"}


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
        return "04.1 Watterson's Theta"
    elif stat == "pi":
        return "04.2 Nucleotide Diversity (Pi)"
    elif stat == "tajima":
        return "04.3 Tajima's D"
