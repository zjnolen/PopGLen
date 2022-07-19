# Genotype likelihood population genomics pipeline

**Under development - mostly works, but not on all systems. And certainly 
needs some documentation... if you'd like to use in this early state, please 
contact me.**

This pipeline performs multiple common population genomic analyses within a 
genotype likelihood framework primarily using ANGSD and related softwares. It 
requires sequencing reads and a reference genome. Right now, it only works 
with libraries with a single run of paired end reads per sample, but this will 
be improved down the line. Some analyses are only possible on UPPMAX's Rackham 
cluster, but will be updated to use singularity as the pipeline matures. As 
the focus is primarily on conservation genomics projects, these will be the 
analyses added first.

## Features

Currently, the pipeline performs the following tasks:

**Reference genome preparation**

- Indexing of reference for subsequent analyses

**Raw read preparation**

- Trimming of paired reads from high quality libraries
- Collapsing of paired reads from fragmented libraries

**Read mapping**

- Mapping prepared raw reads to reference using bwa-mem and clipping of 
  overlapping reads
- Removal of PCR and sequencing duplicates separately for high quality 
  (Picard) and fragmented (DeDup) reads
- Realignment around indels
- Indexing of deduplicated, realigned, mapped reads

**Sample quality control**

- Assess post-mortem DNA damage with DamageProfiler
- Assess mapping quality stats with Qualimap
- Assess endogenous content before duplicate reads are removed (need to 
  re-implement)
- Assess sample duplication and relatedness with ngsrelate

**Data quality filtering**

- Analyses can be set with minimum mapping and base quality thresholds
- Exclusion of entire scaffolds (i.e. sex-linked, low quality) through user 
  config (both list and size based)
- Exclusion of repeat regions from analyses using RepeatModeler/RepeatMasker
- Exclusion of low mappability regions with GenMap
- Exclusion of sites with extreme global depth values (determined separately 
  for the entire dataset, low coverage, and high coverage subsets, then merged)
- Exclusion of sites based on data missingness

**GL based population genomic analyses**

To speed up the pipeline, many of these analyses are done for part of the 
genome at a time, then later merged. This is only done for analyses where 
possible and where the time saved is appreciable. These chunks are made to be 
a user configured size to allow tuning of run-times (i.e. more jobs, shorter 
runtimes vs fewer jobs, longer runtimes).

SAF based analyses are done on variable and non-variable sites passing quality 
filters. This set is the same across all populations in the dataset. Beagle 
based analyses are done on a SNP set that is constant across all populations, 
determined from the output of the Beagle file for the whole dataset. Pruning 
is done on the whole dataset beagle file and the same pruned sites are used 
for all populations.

Additionally, all analyses can be repeated with samples subsampled to a lower 
user configured depth. This helps to ensure results are not simply due to 
variance in depth between groups.

- Linkage pruning where relevant with ngsLD
- PCA with PCAngsd
- Admixture with NGSAdmix
- 1D and 2D Site frequency spectrum production with ANGSD
- Neutrality statistics per population (Watterson's theta, pairwise pi, 
  Tajima's D) in user defined sliding windows with ANGSD
- Estimation of heterozygosity per sample from 1D SFS with ANGSD
- Pairwise Fst between all populations with ANGSD
- Inbreeding coefficients and runs of homozygosity per sample with ngsF-HMM 
  (**NOTE** This is currently only possible for samples that are within a 
  population sampling currently, not for lone samples)

**Reporting of results**

A Snakemake report can be generated after a successful run. The output is 
rather rough for now, but has some improvements in mind:

- Categorization of outputs
- Summary tables of mapping and filtering stats

### Planned

Some additional components to the pipeline are planned, the order below 
roughly corresponding to priority:

- Adjust quality scores of damaged bases with MapDamage (for now, this can 
  still be handled roughly by adding the `-trim INT` or `-noTrans 1` options 
  in the ANGSD parameters)
- Improve portability of the following using containers:
  - PCAngsd
  - ngsF-HMM
  - ngsrelate
- Add options to include genotype calling based approaches for comparison?

## Stages

The pipeline is written to be performed in stages - points at which it can be 
run to and then checked to ensure everything is working properly thus far and 
settings adjusted as needed. **NOTE: This section was written awhile ago - it 
is not up to date.**

### 1) Sequence data and reference processing

This is the first stage and these items are largely processed in parallel. 
This includes the following for each:

**Reference genome**

The reference genome is indexed to prepare it for input into various 
softwares. Additionally, a filter is built for the analyses based on the 
following (all can be turned on or off):

* Scaffold size - Scaffolds shorter than a specified length are excluded from 
  the analyses.
* Autosomes - Sex chromosome related scaffolds (and other selected scaffolds) 
  are removed from the analyses.
* Mappability - Sites with low mappability are removed from the analyses.
* Repeats - Repetitive regions are removed from analyses.
* Extreme depth - Regions of low/high depth are calculated for the entire 
  dataset and separately for low and high coverage samples if specified, 
  removing regions that are extreme in any set.
* Excess heterozygosity - Sites that deviate from HWE and are heterozygous 
  across most of the dataset are removed from analyses.

Most of these will catch similar regions - high repeatability regions will 
have low mappability, excess heterozygosity, and extreme depth. It is 
recommended that if you turn any of these off, to at least leave removal of 
repetitive regions on.

**Sequenceing data**

Reads are quality trimmed and overlapping reads are collapsed per sample. 
Trimmed and collapsed reads are mapped to the reference genome and duplicates 
are removed. Endogenous content is calculated and, if present, low coverage 
samples have post-mortem DNA damage assessed. If requested, base quality 
recalibration can be applied to low coverage samples.

#### What to check after stage 1?

Here, check that the reference filters remove a sensible amount of sites. In 
the case that a filter removed too many or two few sites compared to what 
you'd expect, check the logs to ensure everything ran properly.

For samples, now is a good time to assess sample quality and remove samples 
that are of lower quality or misidentified. If you find that you have a high 
amount of post-mortem DNA damage in historical samples, consider enabling 
DNA damage base quality recalibration or using a filter in downstream analyses 
(`-noTrans 1` or `-trim $bp` as extra options in ANGSD).

Finally, if you have samples of varying depth, you can assess how low 
depth your samples get and determine if you want to downsample higher 
depth samples to match the lower depth ones. This will allow you to 
tease out whether a difference between sample depth is what is leading to 
your results. If you enable this, all analyses will be performed with samples 
at their full depth and downsampled to the specified average depth. This can 
also be added on at any time down the line.

### 2) 'Genotyping'

Genotyping is not exactly the best name for this stage, but this is where 
genotype likelihoods are calculated and either output to a file or output to 
a summary statistic. These are only produced for sites passing reference 
filters with reads and bases that pass quality filters. These aren't manually 
enabled, but rather are made as analyses that need them are enabled. Most of 
these will occur in parallel.

**All sites**

All sites passing filters will have allele frequencies calculated for each 
population and sample using genotype likelihoods and stored in a saf file.

**Varible sites**

Genotype likelihoods for variable sites will be calculated and stored for each 
population and the dataset as a whole in Beagle files.

**Pruned variable sites**

Beagle files for all variable sites will be linkage pruned to produce beagle 
files containing only unlinked positions for analyses that require unlinked 
sites.

**Haplotype calls of variable sites**

Single read sampling approaches allow for a random sampling of an allele at 
each to be used as input data. This produces a plink file containing 
haplotypes for all individuals produced through sampling an allele an allele 
from reads and bases that pass quality filters.

#### What to check after stage 2?

Stage 2 occurs when an analyses needs to be performed, so checking only needs 
to occur when an analysis happens. Once created, summary tables will be made 
that summarize how many sites are included in each of these files - check to 
ensure these make sense, they shouldn't be too much lower than what you have 
from your reference filter.

### 3) Analyses

Analyses can all mostly be run in parallel and are enabled and customized in 
the config file. The available analyses as of now are as follows. Some will 
require additional information to perform.

**PCA with PCAngsd**

**Admixture with NGSAdmix**

**Neutrality statistics and heterozygosity with ANGSD**

* For estimation of neutrality statistics, samples must be separated into 
  populations.

**Inbreeding and Runs of Homozygosity (ROH) with ngsF-HMM**

**Estimation of Effective Migration Surfaces with EEMS**

* Samples must be georeferenced