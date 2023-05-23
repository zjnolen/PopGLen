# Genotype likelihood population genomics pipeline

This workflow aimmed at processing raw sequencing reads and calculating population
genomic statistics within a genotype likelihood framework. As it is focused on GL
methods, it has options to adapt the workflow for processing data with low or variable
coverage and/or contains historical/ancient samples with degraded DNA. It is under
active development, so new features will be added.

## Getting Started

To run this workflow, you'll need paired-end raw sequencing data and a reference genome
to map it to. Currently, the workflow is only compatiblity with datasets containing only
a single library and paired-end sequencing run per sample, and reference genomes must be
uncompressed.

To run the workflow, you will need to be working on a machine with the following:

- Conda (or, preferably, it's drop-in replacement mamba)
- Singularity/Apptainer

Once these softwares are installed, to deploy and configure the workflow, you can follow
the instructions provided in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=zjnolen/angsd-snakemake-pipeline).
You can refer to the Snakemake Documentation for additional information that may be
relevant to your computing environment (running jobs through cluster job queues, setting
default resources).

## Features

Currently, the pipeline performs the following tasks:

### Reference genome preparation

- Indexing of reference for subsequent analyses

### Raw read preparation

- Trimming of paired-end reads from high quality libraries
- Collapsing of paired-end reads from fragmented (aDNA/historical DNA) libraries

### Read mapping

- Mapping prepared raw reads to reference using bwa-mem and clipping of overlapping
  reads
  - **NOTE**: Reads marked as historical (degraded) will only map reads short reads that
    overlap and collapse, to reduce mapping of likely contaminants.
- Removal of PCR and sequencing duplicates separately for fresh (Picard) and fragmented
  (DeDup) DNA reads
- Realignment around indels
- Optional recalibration of base quality scores on degraded DNA bam files with
  [MapDamage2](https://ginolhac.github.io/mapDamage/)
- Indexing of deduplicated, realigned, mapped, and recalibrated reads

### Sample quality control

- Assess post-mortem DNA damage with DamageProfiler
- Assess mapping quality stats with Qualimap
- Assess endogenous content using mapping proportion before duplicate reads are removed

### Data quality filtering

- Analyses can be set with minimum mapping and base quality thresholds
- Exclusion of entire scaffolds (i.e. sex-linked, low quality) through user config (both
  list and contig size based)
- Exclusion of repeat regions from analyses using RepeatModeler/RepeatMasker
- Exclusion of low mappability regions with GenMap
- Exclusion of sites with extreme global depth values (determined separately for the
  entire dataset, and subsets at certain coverage ranges, then merged)
- Exclusion of sites based on data missingness across dataset and/or per population

### GL based population genomic analyses

To speed up the pipeline, many of these analyses are done for part of the genome at a
time, then later merged. This is only done for analyses where possible and where the
time saved is appreciable. These chunks are made to be a user configured size to allow
tuning of run-times (i.e. more jobs, shorter runtimes vs fewer jobs, longer runtimes).

SAF based analyses are done on variable and non-variable sites passing quality filters.
This set is the same across all populations in the dataset and is based on the positions
passing all the requested filters. Beagle (SNP) based analyses are done on a SNP set
that is constant across all populations, determined from the output of the Beagle file
for the whole dataset, and major and minor alleles are inferred from the whole dataset.
When relevant, pruned SNPs are used. Pruning is done on the whole dataset beagle file
and the same pruned sites are used for all populations.

Additionally, all analyses can be repeated with samples subsampled to a lower user
configured depth. This helps to ensure results are not simply due to variance in depth
between groups.

**Analyses:**

- Linkage pruning where relevant with ngsLD
- PCA with PCAngsd
- Admixture with NGSAdmix
- Relatedness using NgsRelate and methods from Waples et al. 2019, *Mol. Ecol.*
- 1D and 2D Site frequency spectrum production with ANGSD
- Neutrality statistics per population (Watterson's theta, pairwise pi, Tajima's D) in
  user defined sliding windows with ANGSD
- Estimation of heterozygosity per sample from 1D SFS with ANGSD
- Pairwise $F_{ST}$ between all populations or individuals in user defined sliding
  windows with ANGSD
- Inbreeding coefficients and runs of homozygosity per sample with ngsF-HMM (**NOTE**
  This is currently only possible for samples that are within a population sampling,
  not for lone samples which will always return an inbreeding coefficient of 0)

### Planned

Some additional components to the pipeline are planned, the order below roughly
corresponding to priority:

- Allow custom site list - either as a supplement to the filters already present or as
  the only filter (by setting all other filters to `false`)
- Add calculation of bootstrapped SFS
- Estimate LD decay with ngsLD
- Manhattan plots in report for sliding window results
- Allow starting with bam files - for those that want to process raw reads in their own
  way before performing analyses
- Add calculation of Dxy
- Add schema for configuration files to improve incorrect format handling and to enable
  defaults

## Workflow directed action graph

Here is a rough graph of the workflow with some rules removed for readability.
It will at some point get a readability improvement when the workflow is
finalized! This is also not really how it looks in the current version, but
it is a close approximation of the flow as of now.

![A directed action graph (DAG) of the main workflow](dag.svg)

## Acknowledgements

The computations required for developing and testing this workflow has been enabled by
resources provided by the National Academic Infrastructure for Supercomputing in Sweden
(NAISS) and the Swedish National Infrastructure for Computing (SNIC) at UPPMAX partially
funded by the Swedish Research Council through grant agreements no. 2022-06725 and no.
2018-05973.
