# Genotype likelihood population genomics pipeline

This workflow is aimed at processing sequencing data and calculating population
genomic statistics within a genotype likelihood framework. As a primary use
case of genotype likelihood based methods is analysis of samples sequenced to
low depth or from degraded samples, processing can optionally be adapted to
account for DNA damage. It is aimed at single and multi-population analyses of
samples mapped to the same reference, so is appropriate for datasets containing
individuals from a single or multiple populations and allows for examining
population structure, genetic diversity, genetic differentiation, allele
frequencies, linkage disequilibrium, and more.

The workflow is designed with two entry points in mind. Users with raw
sequencing data in FASTQ format stored locally or as an NCBI/ENA run accession
can perform raw sequencing data alignment and processing, followed up by the
population genomic analyses. Alternatively, users who have already mapped and
processed their reads into BAM files can use these to start directly at the
population genomic analyses (for instance, if you bring BAM files from
GenErode, nf-core/eager, or your own custom processing).

## Features

### Sequence data processing

If starting with raw sequencing data in FASTQ format, the pipeline will handle
mapping and processing of the reads, with options speific for historical DNA.
These steps are available if providing paths to local FASTQ files or SRA
accessions from e.g. NCBI/ENA.

- Trimming and collapsing of overlapping paired end reads with fastp
- Mapping of reads using bwa mem
- Estimation of read mapping rates (for endogenous content calculation)
- Removal of duplicates using Picard for paired reads and dedup for collapsed
  overlapping reads
- Realignment around indels using GATK
- Optional base quality recalibration using MapDamage2 in historical samples to
  correct post-mortem damage
- Clipping of overlapping mapped paired end reads using bamutil
- Quality control information from fastp, Qualimap, DamageProfiler, and
  MapDamage2 (Qualimap is also available for users starting with BAM files)

### Population Genomics

The primary goal of this pipeline is to perform analyses in ANGSD and related
softwares in a repeatable and automated way. These analyses are available both
when you start with raw sequencing data or with processed BAM files.

- Estimation of linkage disequilibrium across genome and LD decay using ngsLD
- Linkage pruning where relevant with ngsLD
- PCA with PCAngsd
- Admixture with NGSAdmix
- Relatedness using NgsRelate and IBSrelate
- 1D and 2D Site frequency spectrum production with ANGSD
- Neutrality statistics per population (Watterson's theta, pairwise nucleotide
  diversity, Tajima's D) in user defined sliding windows with ANGSD
- Estimation of heterozygosity per sample from 1D SFS with ANGSD
- Pairwise $F_{ST}$ between all populations or individuals in user defined
  sliding windows with ANGSD
- Inbreeding coefficients and runs of homozygosity per sample with ngsF-HMM
  (**NOTE** This is currently only possible for samples that are within a
  population sample, not for lone samples which will always return an
  inbreeding coefficient of 0)
- Identity by state (IBS) matrix between all samples using ANGSD

Additionally, several data filtering options are available:

- Identification (can also skip identification if repeat bed/gff is supplied)
  and removal of repetitive regions
- Removal of regions with low mappability for fragments of a specified size
- Removal of regions with extreme high or low depth
- Removal of regions with a certain amount of missing data
- Multiple filter sets from user provided BED files that can be intersected
  with other enabled filters (for instance, performing analyses on neutral
  sites and genic regions separately)

All the above analyses can also be performed with sample depth subsampled to
a uniform level to account for differences in depth between samples.

## Getting Started

The workflow will require data in either FASTQ or BAM format, as well as a
single (uncompressed) reference genome that the reads will be or have been
mapped to. Analyses are intended to be between samples and populations of
samples mapped to the same reference.

To run the workflow, you will need to be working on a machine with the
following:

- Conda (or, preferably, it's drop-in replacement mamba)
- Singularity/Apptainer

Once these softwares are installed, to deploy and configure the workflow, you
can follow the instructions provided in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=zjnolen/angsd-snakemake-pipeline).
You can refer to the Snakemake Documentation for additional information that
may be relevant to your computing environment (running jobs through cluster job
queues, setting default resources).

### Notes on inputs

#### Reference Genomes

Reference genomes should be uncompressed, and contig names should be clear and
concise. Currently, there are some issues parsing contig names with
underscores, so please change these in your reference before running the
pipeline. Alphanumeric characters, as well as `.` in contig names have been
tested to work so far, other symbols have not been tested.

Potentially the ability to use bgzipped genomes will be added, I just need to
check that it works with all underlying tools. Currently, it will for sure not
work, as calculating chunks is hard-coded to work on an uncompressed genome.

#### BAM Files

BAM files will receive no processing, aside from optionally clipping
overlapping read pairs that are otherwise double counted in ANGSD. Ensure any
processing (duplicate removal, damage correction, realignment) you wish to
include has already been performed. Ideally, it would be good to clip the
overlapping reads before running the workflow if you are supplying your own
bams, and turn off clipping of user provided bams in the pipeline config. This
prevents doubling the storage usage, as otherwise all your bams get duplicated,
but with clipped reads. Doing this beforehand with [bamutil](https://genome.sph.umich.edu/wiki/BamUtil)
is straightforward: `bam clipOverlap --in {input.bam} --out {output.bam}`.
After this, you can remove the original bams to save storage.

This pipeline is written with reusing of samples across datasets in mind. For
samples starting at fastq, this should be done seamlessly by reusing samples
in different datasets (as set in the config file) processed in the same working
directory. This, in most cases, will also be true for bam input. However, if
you are planning to process datasets where different bam files will be used
as input for a sample given the same sample name and reference name in both
dataset configs (for instance, to run the same samples with and without
MapDamage rescaling), the second dataset will overwrite the bams of the first,
and Snakemake will suggest rerunning the first dataset. In this case, it is
best to have different sample names for the different input bams, or to run the
two datasets in different working directories.

Some QC will not be available for users starting at BAM files. No read
processing QC can be produced and should be done beforehand. While mapping
percentages are calculated, these may not entirely represent the truth, as they
may not account for anything already fully removed from the bam. In this case,
they also can't be separated into categories of collapsed and uncollapsed
reads, and instead are simply reported as the total percentage mapping only.

### Running on a cluster

Development was done on UPPMAX's Rackham cluster, and a simple profile is
included in the [`rackham`](rackham) folder to simplify running this workflow
through SLURM there. For running on other SLURM based cluster configs, this
file should largely work with a few minor modifications of the defaults.
Primarily, this means ensuring that the resources make sense for your system,
i.e. changing the default partition and account. Not changing the default
memory may result in over-reserving memory in some cases, but a quick fix would
be to change `6400` to whatever the default memory reserved per cpu is on your
HPC (though then you might need to up the threads requested fo in some rules).
See [Snakemake's cluster support documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
for information on how to adapt the profile for your HPC environment.

#### Resources (and what to do if the workflow is over/underbooking resources)

As Rackham ties memory reservations to cpu reservations, the resource
allocation for rules is mostly done through threads currently. In the future
thread and memory resources will be more explicitly defined. For now, if you
find that the workflow is over/underbooking a given resource, you can adjust
the resource reservations in your run profile. See the commented section in
the [`rackham/config.yaml`](rackham/config.yaml) config to see an example of
this. Right now, memory is only ever defined through threads, so you may need
to lower the threads and add a memory resource to some rules using this method
in order to optimize them for your system.

## Workflow directed action graph

Below is a graph of the workflow with most of the analyses enabled. This is
generated directly by Snakemake, though I have removed the `all`, `popfile`,
and `link_ref` rules to improve readability, as they are not needed to
understand the analysis flow. A more condensed, readable diagram will be added
shortly.

In general, there a few stages that can be seen grouping here. A mapping stage
at the top, ending with `symlink_bams`, followed by the creation of a set of
filters for the dataset that ends with `combine_beds`. Afterwards, population
genetic analyses are performed in primarily two pathways - allele frequency
based results (starting with `angsd_doSaf_sample/pop`) and SNP based results
(starting with `angsd_doGlf2` (beagle)).

![A directed action graph (DAG) of the main workflow](images/rulegraph.svg)

## Acknowledgements

The computations required for developing and testing this workflow has been
enabled by resources provided by the National Academic Infrastructure for
Supercomputing in Sweden (NAISS) and the Swedish National Infrastructure for
Computing (SNIC) at UPPMAX partially funded by the Swedish Research Council
through grant agreements no. 2022-06725 and no. 2018-05973.
