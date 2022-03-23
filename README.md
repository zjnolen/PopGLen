# Genotype likelihood population genomics pipeline

This pipeline focuses on performing population genomic analyses within a 
genotype likelihood framework primarily using ANGSD and related softwares. It 
requires sequencing reads and a reference genome. Right now, it only works 
with libraries with a single set of paired end reads per sample, but this will 
be improved down the line. Some analyses are only possible on UPPMAX's Rackham 
cluster, but will be updated to conda and singularity as the pipeline matures. 
As the focus is primarily on conservation genomics projects, these will be the 
analyses added first.

### Currently implemented

#### Mapping

- Fastq processing with fastp - trims low quality bases, poly-g tails, and 
  collapses overlapping read pairs
- Mapping to reference with bwa-mem - maps both collapsed and uncollapsed 
  paired reads and singletons
- Endogenous content estimation with samtools flagstat
- Duplicate removal with Picard
- Mapping quality control with Qualimap

#### Analyses

- PCA with PCAngsd (only on Rackham for now)
- NGSadmix (needs publishable wrapper, GitHub version won't work until this is 
  done)
- Estimation of thetas (pi, Watterson's, Tajima's D) with ANGSD
- Estimation of per individual heterozygosity with ANGSD (estimates SFS, but
  full estimation of values not ready)
- Estimation of inbreeding coefficients and runs of homozygosity with ngsF-HMM 
  (only on Rackham for now)

### Planned updates

- GONE (not GL based)

### Current organization

Things should mostly be runnable now, working on documentation. Master branch 
only has fully functioning commits for currently implemented items, develop is 
where new stuff is added incrementally.