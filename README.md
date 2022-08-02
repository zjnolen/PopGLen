# Genotype likelihood population genomics pipeline

**Under development - mostly works, but not on all systems. And certainly 
needs some documentation... if you'd like to use in this early state, give it 
a try, and feel free to send me questions.**

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

- Trimming of paired-end reads from high quality libraries
- Collapsing of paired-end reads from fragmented libraries

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
- ~~Assess endogenous content before duplicate reads are removed~~ Needs to be 
  re-added
- Assess sample duplication and relatedness with NgsRelate

**Data quality filtering**

- Analyses can be set with minimum mapping and base quality thresholds
- Exclusion of entire scaffolds (i.e. sex-linked, low quality) through user 
  config (both list and size based)
- Exclusion of repeat regions from analyses using RepeatModeler/RepeatMasker
- Exclusion of low mappability regions with GenMap
- Exclusion of sites with extreme global depth values (determined separately 
  for the entire dataset, low coverage, and high coverage subsets, then merged)
- Exclusion of sites based on data missingness across dataset

**GL based population genomic analyses**

To speed up the pipeline, many of these analyses are done for part of the 
genome at a time, then later merged. This is only done for analyses where 
possible and where the time saved is appreciable. These chunks are made to be 
a user configured size to allow tuning of run-times (i.e. more jobs, shorter 
runtimes vs fewer jobs, longer runtimes).

SAF based analyses are done on variable and non-variable sites passing quality 
filters. This set is the same across all populations in the dataset and is 
based on the positions passing all the requested filters. Beagle (SNP) based analyses are done on a SNP set that is constant across all populations, 
determined from the output of the Beagle file for the whole dataset. When 
relevant, pruned SNPs are used. Pruning is done on the whole dataset beagle 
file and the same pruned sites are used for all populations.

Additionally, all analyses can be repeated with samples subsampled to a lower 
user configured depth. This helps to ensure results are not simply due to 
variance in depth between groups.

**Analyses:**
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
  population sampling, not for lone samples which will always return an 
  inbreeding coefficient of 0)

**Reporting of results**

A Snakemake report can be generated after a successful run. The output is 
rather rough for now, but has some improvements in mind.

### Planned

Some additional components to the pipeline are planned, the order below 
roughly corresponding to priority:

**Data analysis**
- Allow custom site list - either as a supplement to the filters already 
  present or as the only filter (by setting all other filters to `false`)
- Allow starting with bam files - for those that want to process raw reads 
  in their own way before performing analyses
- Adjust quality scores of damaged bases with MapDamage (for now, this can 
  still be handled roughly by adding the `-trim INT` or `-noTrans 1` options 
  in the ANGSD parameters)
- Add calculation of Dxy (either Peñalba's or Marques's method)
- Improve portability of the following using containers:
  - PCAngsd
  - ngsF-HMM
  - ngsrelate
- TreeMix using Maf file outputs?
- Add options to include genotype calling based approaches for comparison?

**Reporting**

- Categorization of outputs
- Summary tables of mapping and filtering stats

## Stages

It is most often best to perform the pipeline in stages, check that things 
are working up to that point, then add more in. This makes for simpler 
troubleshooting and also prevents needless computation time on downstream 
analyses in the event of a mistake earlier in the pipeline. How much of the 
pipeline you run is largely decided by the `analyses` section of the config 
file. After each successful run, you can create a report with the results thus 
far, helping you to decide if you wish to add more. I suggest starting with everything as `false` and adding analyses in the following order:

### 1) Sequence data and reference processing, quality control

This is the first stage, where your reference and raw data are processed and 
filtered dataset developed. To initiate this stage, set all the filters you 
intend to use under the `# filtering` section of `analyses` in the config 
file to `true`. Set any quality control outputs you'd like to output to `true` 
as well. If you've included `relatedness`, this is enough for the first run. 
If not, set one downstream analysis that uses all the data to `true`,  to 
ensure all your files and filters get processed. `pca_pcangsd` is a quick and 
simple one to list as `true` and ensure everything runs. It also is a simple 
result to look at and see if anything might be off.

#### Checks afterwards

After this initial run, check out the filtering summary in the 
`results/{dataset}/genotyping/filters/beds/{dataset}_filts.sum` file. This 
shows how many sites pass each filter and what percent of the genome that 
covers. If this seems off, take a look at the settings for your filtering and adjust.

Next, check out the quality reports of individual samples. Should any related 
individuals be removed based on the NgsRelate output? Individuals with 
exceptionally low mapping quality or coverage based on the Qualimap outputs? 
Do any individuals have DNA damage in need of correction based on the 
DamageProfiler outputs? If you add these samples to `exclude_ind` in the 
config and re-run the pipeline, Snakemake will re-run rules that included 
these individuals without them.

### 2) Analyses

If everything looks good after this initial stage, you can set any remaining 
analyses you would like to `true` in the config. At this stage you can also 
determine a downsampling coverage to include in the `downsample_cov` option if 
you wish to perform all analyses a second time with samples all downsampled to 
the same coverage. You can use the Qualimap coverages to determine what level 
you want to set this too. Any analysis can be added on later, so feel free to 
add as few or as many as you prefer at this point.

#### Checks afterwards

This should be your final results, so ensure your outputs in the report for 
each analysis make sense and adjust settings for relevant analyses. Adjusting 
parameters will be noticed by snakemake, and simply re-running the command 
will result in re-running of relevant portions of the pipeline. This is true 
for the addition and removal of samples from the sample list as well.