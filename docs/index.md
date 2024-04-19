# Welcome to the documentation for ANGSD Snakemake Pipeline

ANGSD Snakemake Pipeline is aimed at enabling users to run population genomic
analyses on their data within a genotype likelihood framework in an automated
and reproducible fashion. Genotype likelihood based analyses avoid genotype
calling, instead performing analyses on the likelihoods of each possible
genotype, incorporating uncertainty about the true genotype into the analysis.
This makes them especially suited for datasets with low coverage or that vary
in coverage.

This pipeline was developed in large part to make my own analyses easier. I
work with many species being mapped to their own references within the same
project. I developed this pipeline so that I could ensure standardized
processing for datasets within the same project and to automate the many steps
that go into performing these analyses. As it needed to fit many datasets, it
is generalizable and customizable through a single configuration file and uses
a common workflow utilized by ANGSD users, so it is available for others to
use, should it suit their needs.

??? question "Questions? Feature requests? Just ask!"
    I'm glad to answer questions on the [GitHub Issues](https://github.com/zjnolen/angsd-snakemake-pipeline/issues)
    page for the project, as well as take suggestions for features or
    improvements!

## Pipeline Summary

The pipeline aims to follow the general path many users will use when working
with ANGSD and other GL based tools. Raw equencing data is processed into BAM
files (with optional configuration for historical degraded samples) or BAM
files are provided directly. From there several quality control reports are
generated to help determine what samples are included. The pipeline then builds
a 'sites' file to perform analyses with. This sites file is made from several
user-configured filters, intersecting all and outputing a list of sites for the
analyses to be performed on across all samples. This can also be extended by
user-provided filter lists (e.g. to limit to neutral sites, genic regions,
etc.).

After samples have been processed, quality control reports produced, and the
sites file has been produced, the pipeline can continue to the analyses.

- Linkage disequilibrium estimation, LD decay, LD pruning (ngsLD)
- PCA (PCAngsd)
- Admixture (NGSAdmix)
- Inbreeding/Runs of Homozygosity (ngsF-HMM)
- Relatedness (NGSRelate, IBSrelate)
- Identity by state matrix (ANGSD)
- Site frequency spectrum (ANGSD)
- Watterson's estimator ($θ_w$), Nucleotide diversity ($π$), Tajima's $D$
  (ANGSD)
- Individual heterozygosity with bootstrapped confidence intervals (ANGSD)
- Pairwise $F_{ST}$ (ANGSD)

These all can be enabled and processed independently, and the pipeline will
generate genotype likelihood input files using ANGSD and share them across
analyses as appropriate, deleting temporary intermediate files when they are no
longer needed.

At any point after a successful completion of a portion of the pipeline, a
report can be made that contains tables and figures summarizing the results
for the currently enabled parts of the pipeline.

If you're interested in using this, head to the
[Getting Started](getting-started.md) page!
