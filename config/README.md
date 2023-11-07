# General Settings

In the [`config`](./) directory, you will find template configuration files
for this pipeline. For your run, you'll need to edit `samples.tsv`,
`units.tsv`, and `config.yaml`.

## Samples list

All your samples should be listed in `samples.tsv`. Samples preceded with a `#`
will not be included, this can be useful if you want to exclude a sample after
quality checking.

Each sample must have four columns filled in the sample sheet. The columns are
tab separated:

- `sample` - The name of the sample
- `population` - The populations you will group your samples into. These are
  the groups that population level stats are calculated on.
- `time` - This should be either `modern` or `historical`, the only thing this
  will affect is whether or not your bam files will be corrected for
  post-mortem damage or not.
- `depth` - This is only used for depth filtering. Extreme levels (both high
  and low) will be calculated for each group you put here as well as the
  dataset as a whole, and all will be filtered out for all analyses. Any string
  can be used for this. If all samples are sequenced to roughly similar depths,
  all can have the same value. If some are low coverage, and some higher,
  simply using 'low' and 'high' on the corresponding samples is sufficient.

The values in the sample list will end up in filenames, so ensure that they are
only characters permitted in filenames on your system.

## Units list

All your raw data will be pointed to in `units.tsv`.

Each sample must have five columns filled in the units sheet. The columns are
tab separated:

- `sample` - The sample name, same as in `samples.tsv`.
- `unit` - This is used to fill out the `ID` read group in the bam file. It
  must be unique to each read group, so the same sample shouldn't have the
  same unit for more than one sequencing run. A good format might be
  `sequencer_barcode.lane`. Optical duplicates will be removed within units.
- `lib` - This is used to fill out the `LB` read group. This should be a unique
  identifier for each sample library. Sequencing runs from the same library,
  but different runs, will have the same value in `lib`, but different in
  `unit`.
- `platform` - This is used to fill out the `PL` read group. Put what you'd
  want there. Usually `ILLUMINA` for Illumina platforms.
- `fq1` and `fq2` - The absolute or relative paths from the working directory
  to the raw fastq files for the sample. Currently the pipeline only supports
  paired-end sequencing, single end may be added down the line.

## Configuration file

`config.yaml` contains the configuration for the workflow, this is where you
will put what analyses, filters, and options you want. Below I describe the
configuration options. The [`config.yaml`](config.yaml) in this repository
serves as a template, but includes some 'default' parameters that may be good
starting points for some users. If `--configfile` is not specified in the
snakemake command, the workflow will default to `config/config.yaml`.

### Configuration options

#### Dataset Configuration

Required configuration of the 'dataset'.

- `samples:` An absolute or relative path from the working directory to the
  `samples.tsv` file.
- `units:` An absolute or relative path from the working directory to the
  `units.tsv` file.
- `dataset:` A name for this dataset run - essentially, an identifier for a
  batch of samples to be analysed together with the same configuration.

Here, dataset means a set of samples and configurations that the workflow will
be run with. Each dataset should have its own `samples.tsv` and `config.yaml`,
but the same `units.tsv` can be used for multiple if you prefer. Essentially,
what the dataset identifier does is keeps your outputs organized into projects,
so that the same BAM files can be used in multiple datasets without having to
be remade.

So, say you have `dataset1_samples.tsv` and `dataset2_samples.tsv`, with
corresponding `dataset1_config.tsv` and `dataset2_config.yaml`. The sample
files contain different samples, though some are shared between the datasets.
The workflow for dataset1 can be run, and then dataset2 can be run. When
dataset2 runs, it map new samples, but won't re-map samples processed in
dataset1. Each will perform downstream analyses independently with their sample
set and configuration files, storing these results in dataset specific folders.

#### Reference Configuration

Required configuration of the reference.

- `chunk_size:` A size in bp (integer). Your reference will be analyzed in
  'chunks' of contigs of this size to parallelize processing. This size should
  be larger than the largest contig in your genome. A larger number means fewer
  jobs that run longer. A smaller number means more jobs that run shorter. The
  best fit will depend on the reference and the compute resources you have
  available. Leaving this blank will not divide the reference up into chunks
  (but this isn't optimized yet, so it will do a couple unnecessary steps).

- `reference:`
  - `name:` A name for your reference genome, will go in the file names.
  - `fasta:` A path to the reference fasta file (currently only supports
    uncompressed fasta files)
  - `mito:` Mitochondrial contig name(s), will be removed from analysis. Should
    be listed within brackets []
  - `sex-linked:` Sex-linked contig name(s), will be removed from analysis.
    Should be listed within brackets []
  - `exclude:` Additional contig name(s) to exclude from analysis. Should be
    listed within brackets []
  - `min_size:` A size in bp (integer). All contigs below this size will be
    excluded from analysis.

Reference genomes should be uncompressed, and contig names should be clear and
concise. Currently, there are some issues parsing contig names with
underscores, so please change these in your reference before running the
pipeline. Alphanumeric characters, as well as `.` in contig names have been
tested to work so far, other symbols have not been tested.

Potentially the ability to use bgzipped genomes will be added, I just need to
check that it works with all underlying tools. Currently, it will for sure not
work, as calculating chunks is hard-coded to work on an uncompressed genome.

#### Sample Set Configuration

This will exclude individuals from analysis that are listed in the sample list.
This may be useful if you run the workflow and find a poor quality sample, and
want to re-run without it. Or if you have relatives in the dataset and you want
to exclude them where necessary.

- `exclude_ind:` Sample name(s) that will be excluded from the workflow. Should
  be a list in [].
- `excl_pca-admix:` Sample name(s) that will be excluded *only* from PCA and
  Admixture analyses. Useful for close relatives that violate the assumptions
  of these analyses, but that you want in others. Should be a list in [].

#### Analysis Selection

Here, you will define which analyses you will perform. It is useful to start
with only a few, and add more in subsequent workflow runs, just to ensure you
catch errors before you use compute time running all analyses. Most are set
with (`true`/`false`) or a value, described below. Modifications to the
settings for each analysis are set in the next section.

- `populations:` A list of populations found in your sample list to limit
  population analyses to. Might be useful if you want to perform individual
  analyses on some samples but not include them in any population level
  analyses

- `analyses:`
  - `genmap:` Filter out sites with low mappability estimated by Genmap
  (`true`/`false`)
  - `repeatmasker:` (NOTE: Only one of the three options should be filled/true)
    - `local_lib:` Filter repeats by masking with an already made library you
      have locally. Should be file path.
    - `dfam_lib:` Filter repeats using a library available from dfam. Should be
      a taxonomic group name.
    - `build_lib:` Use RepeatModeler to build a library of repeats from the
      reference itself, then filter them from analysis (`true`/`false`).
  - `extreme_depth:` Filter out sites with extremely high or low global
    sequencing depth (`[lower, upper]`). The value of `lower` (float) will be
    multiplied by the median global depth to create a lower global depth
    threshold, and `upper` will do the same to creat an upper threshold. This
    is done for all samples, as well as separately for depth groupings defined
    in samples file.
  - `dataset_missing_data:` A floating point value between 0 and 1. Sites with
    data for fewer than this proportion of individuals across the whole dataset
    will be filtered out.
  - `population_missing_data:` A floating point value between 0 and 1. Sites
    with data for fewer than this proportion of individuals in any population
    will be filtered out for all populations.
  - `qualimap:` Perform Qualimap bamqc on bam files for general quality stats
    (`true`/`false`)
  - `damageprofiler:` Estimate post-mortem DNA damage on historical samples
    with Damageprofiler (`true`/`false`) NOTE: This just adds the addition of
    Damageprofiler to the already default output of MapDamage. DNA damage will
    always be estimated and rescaled by MapDamage for samples marked as
    'historical'
  - `estimate_ld:` Estimate pairwise linkage disquilibrium between sites with
    ngsLD for each popualation and the whole dataset. Note, only set this if
    you want to generate the LD estimates for use in downstream analyses
    outside this workflow. Other analyses within this workflow that require LD
    estimates (LD decay/pruning) will function properly regardless of the
    setting here. (`true`/`false`)
  - `ld_decay:` Use ngsLD to plot LD decay curves for each population and for
    the dataset as a whole (`true`/`false`)
  - `pca_pcangsd:` Perform Principal Component Analysis with PCAngsd
    (`true`/`false`)
  - `admix_ngsadmix:` Perform admixture analysis with NGSadmix (`true`/`false`)
  - `relatedness:` Can be performed multiple ways, set any combination of the
    three options below. Note, that I've mostly incorporated these with the
    R0/R1/KING kinship methods in Waples et al. 2019, *Mol. Ecol.* in mind.
    These methods differ slightly from how they implement this method, and will
    give slightly more/less accurate estimates of kinship depending on your
    reference's relationship to your samples. `ibsrelate_ibs` uses the
    probabilities of all possible genotypes, so should be the most accurate
    regardless, but can use a lot of memory and take a long time with many
    samples. `ibsrelate_sfs` is a bit more efficient, as it does things in a
    pairwise fashion in parallel, but may be biased if the segregating alleles
    in your populations are not represented in the reference. `ngsrelate` uses
    several methods, one of which is similar to `ibsrelate_sfs`, but may be
    less accurate due to incorporating in less data. In my experience,
    NGSrelate is suitable to identify up to third degree relatives in the
    dataset, but only if the exact relationship can be somewhat uncertain (i.e.
    you don't need to know the difference between, say, parent/offspring and
    full sibling pairs, or between second degree and third degree relatives).
    IBSrelate_sfs can get you greater accuracy, but may erroneously inflate
    kinship if your datset has many alleles not represented in your reference.
    If you notice, for instance, a large number of third degree relatives
    (KING ~0.03 - 0.07) in your dataset that is not expected, it may be worth
    trying the IBS based method (`ibsrelate_ibs`).
    - `ngsrelate:` Co-estimate inbreeding and pairwise relatedness with
      NGSrelate (`true`/`false`)
    - `ibsrelate_ibs:` Estimate pairwise relatedness with the IBS based method
      from Waples et al. 2019, *Mol. Ecol.*. This can use a lot of memory, as
      it has genotype likelihoods for all sites from all samples loaded into
      memory, so it is done per 'chunk', which still takes a lot of time and
      memory. (`true`/`false`)
    - `ibsrelate_sfs:` Estimate pairwise relatedness with the SFS based method
      from Waples et al. 2019, *Mol. Ecol.*. Enabling this can greatly increase
      the time needed to build the workflow DAG if you have many samples. As a
      form of this method is implemented in NGSrelate, it may be more
      efficient to only enable that. (`true`/`false`)
  - `thetas_angsd:` Estimate pi, theta, and Tajima's D for each population in
    windows across the genome using ANGSD (`true`/`false`)
  - `heterozygosity_angsd:` Estimate individual genome-wide heterozygosity
    using ANGSD (`true`/`false`)
  - `fst_angsd:` Estimate pairwise $F_{ST}$ using ANGSD. Set one or both of the
    below options. Estimates both globally and in windows across the genome.
    - `populations:` Pairwise $F_{ST}$ is calculated between all possible
      population pairs (`true`/`false`)
    - `individuals:` Pairwise $F_{ST}$ is calculated between all possible
      population pairs. NOTE: This can be really intensive on the DAG building
      process, so I don't recommend enabling unless you're certain you want
      this (`true`/`false`)
  - `inbreeding_ngsf-hmm:` Estimates inbreeding coefficients and runs of
    homozygosity using ngsF-HMM. Output is converted into an inbreeding measure
    $F_ROH$, which describes the proportion of the genome in runs of
    homozygosity over a certain length. (`true`/`false`)
  - `ibs_matrix:` Estimate pairwise identity by state distance between all
    samples using ANGSD. (`true`/`false`)

#### Filter Sets

By default, this workflow will perform all analyses requested in the above
section on all sites that pass the filters set in the above section. These
outputs will contain `allsites-filts` in the filename and in the report.
However, many times, it is useful to perform an analysis on different subsets
of sites, for instance, to compare results for genic vs. intergenic regions,
neutral sites, exons vs. introns, etc. Here, users can set an arbitrary number
of additional filters using BED files. For each BED file supplied, the contents
will be intersected with the sites passing the filters set in the above
section, and all analyses will be performed additionally using those sites.

For instance, given a BED file containing putatively neutral sites, one could
set the following:

```yaml
filter_beds:
  neutral-sites: "resources/neutral_sites.bed"
```

In this case, for each requested analysis, in addition to the `allsites-filts`
output, a `neutral-filts` (named after the key assigned to the BED file in
`config.yaml`) output will also be generated, containing the results for sites
within the specified BED file that passed any set filters.

More than one BED file can be set, up to an arbitrary number:

```yaml
filter_beds:
  neutral: "resources/neutral_sites.bed"
  intergenic: "resources/intergenic_sites.bed"
  introns: "resources/introns.bed"
```

It may also sometimes be desireable to skip analyses on `allsites-filts`, say
if you are trying to only generate diversity estimates or generate SFS for a
set of neutral sites you supply. To skip running any analyses for
`allsites-filts` and only perform them for the BED files you supply, you can
set `only_filter_beds: true` in the config file. This may also be useful in the
event you have a set of already filtered sites, and want to run the workflow on
those, ignoring any of the built in filter options by setting them to `false`.

#### Software Configuration

These are software specific settings that can be user configured in the
workflow. If you are missing a configurable setting you need, open up an issue
or a pull request and I'll gladly put it in.

- `mapQ:` Phred-scaled mapping quality filter. Reads below this threshold will
  be filtered out. (integer)
- `baseQ:` Phred-scaled base quality filter. Reads below this threshold will be
  filtered out. (integer)

- `params:`
  - `genmap:` Parameters for mappability analysis, see [GenMap's documentation](https://github.com/cpockrandt/genmap/)
    for more details.
    - `K:`
    - `E:`
    - `map_thresh:` A threshold mappability score. Sites with a mappability
      score below this threshold are filtered out if GenMap is enabled.
      (integer/float, 0-1)
  - `fastp:`
    - `extra:` Additional options to pass to fastp trimming. (string)
  - `picard:`
    - `MarkDuplicates:` Additional options to pass to Picard MarkDuplicates.
      `--REMOVE_DUPLICATES true` is recommended. (string)
  - `angsd:` General options in ANGSD, relevant doc pages are linked
    - `gl_model:` Genotype likelihood model to use in calculation
      (`-GL` option in ANGSD, [docs](http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods))
    - `maxdepth:` When calculating individual depth, sites with depth higher
      than this will be binned to this value. Should be fine for most to leave
      at `1000`. (integer, [docs](http://www.popgen.dk/angsd/index.php/Depth))
    - `extra:` Additional options to pass to ANGSD during genotype likelihood
      calculation. This is primarily useful for adding BAM input filters. Note
      that `--remove_bads` and `-only_proper_pairs` are enabled by default, so
      they only need to be included if you want to turn them off. I've also
      found that for some datasets, `-C 50` and `-baq 1` can create a strong
      relationship between sample depth and detected diversity, effectively
      removing the benefits of ANGSD for low/variable depth data. I recommend
      that these aren't included unless you know you need them, and even then,
      I'd recommend plotting `heterozygosity ~ sample depth` to ensure there is
      not any relationship. Since the workflow uses bwa to map, `-uniqueOnly 1`
      doesn't do anything if your minimum mapping quality is > 0. Don't put
      mapping and base quality thresholds here either, it will use the ones
      defined above automatically. Although historical samples will have DNA
      damaged assessed and to some extent, corrected, it may be useful to put
      `-noTrans 1` or `-trim INT` here if you're interested in stricter filters
      for degraded DNA. (string, [docs](http://www.popgen.dk/angsd/index.php/Input#BAM.2FCRAM))
    - `snp_pval:` The p-value to use for calling SNPs (float, [docs](http://www.popgen.dk/angsd/index.php/SNP_calling))
    - `min_maf:` The minimum minor allele frequency required to call a SNP.
      (float, [docs](http://www.popgen.dk/angsd/index.php/Allele_Frequencies))
  - `ngsld:` Settings for ngsLD ([docs](https://github.com/fgvieira/ngsLD))
    - `max_kb_dist_est-ld:` For the LD estimates generated when setting
      `estimate_ld: true` above, set the maximum distance between sites in kb
      that LD will be estimated for (`--max_kb_dist` in ngsLD, integer)
    - `max_kb_dist_decay:` The same as `max_kb_dist_est-ld:`, but used when
      estimating LD decay when setting `ld_decay: true` above (integer)
    - `max_kb_dist_pruning:` The same as `max_kb_dist_est-ld:`, but used when
      linkage pruning SNPs as inputs for PCA, Admix, and Inbreeding analyses.
      Any positions above this distance will be assumed to be in linkage
      equilibrium during the pruning process (integer)
    - `rnd_sample_est-ld:` For the LD estimates generated when setting
      `estimate_ld: true` above, randomly sample this proportion of pairwise
      linkage estimates rather than estimating all (`--rnd_sample` in ngsLD,
      float)
    - `rnd_sample_decay:` The same as `rnd_sample_est-ld:`, but used when
      estimating LD decay when setting `ld_decay: true` above (float)
    - `fit_LDdecay_extra:` Additional plotting arguments to pass to
      `fit_LDdecay.R` when estimating LD decay (string)
    - `fit_LDdecay_n_correction:` When estimating LD decay, should the sample
      size corrected r^2 model be used? (`true`/`false`, `true` is the
      equivalent of passing a sample size to `fit_LDdecay.R` in ngsLD using
      `--n_ind`)
    - `pruning_min-weight:` The minimum r^2 to assume two positions are in
      linkage disequilibrium when pruning (float)
  - `realSFS:` Settings for realSFS
    - `fold:` Whether or not to fold the produced SFS (0 or 1, [docs](http://www.popgen.dk/angsd/index.php/SFS_Estimation))
      **NOTE:** I have not implemented the use of an ancestral reference into
      this workflow, so this should always be set to 1 until I implement this.
    - `sfsboot:` Doesn't work now, but when it does it will produce this many
      bootstrapped SFS per population and population pair (integer)
  - `fst:` Settings for $F_{ST}$ calculation in ANGSD
    - `whichFst:` Determines which $F_{ST}$ estimator is used by ANGSD. With 0
      being the default Reynolds 1983 and 1 being the Bhatia 2013 estimator.
      The latter is preferable for small or uneven sample sizes
      (0 or 1, [docs](http://www.popgen.dk/angsd/index.php/Fst))
    - `win_size:` Window size in bp for sliding window analysis (integer)
    - `win_step:` Window step size in bp for sliding window analysis (integer)
  - `thetas:` Settings for pi, theta, and Tajima's D estimation
    - `win_size:` Window size in bp for sliding window analysis (integer)
    - `win_step:` Window step size in bp for sliding window analysis (integer)
  - `ngsadmix:` Settings for admixture analysis with NGSadmix. This analysis is
    performed for a set of K groupings, and each K has several replicates
    performed. Replicates will continue until a set of N highest likelihood
    replicates converge, or the number of replicates reaches an upper threshold
    set here. Defaults for `reps`, `minreps`, `thresh`, and `conv` can be left
    as default for most.
    - `kvalues:` A list of values of K to fit the data to (list of integers)
    - `reps:` The maximum number of replicates to perform per K (integer)
    - `minreps:` The minimum number of replicates to perform, even if
      replicates have converged (integer)
    - `thresh:` The convergence threshold - the top replicates must all be
      within this value of log-likelihood units to consider the run converged
      (integer)
    - `conv:` The number of top replicates to include in convergence
      assessment. (integer)
    - `extra:` Additional arguments to pass to NGSadmix (for instance,
      increasing `-maxiter`). (string, [docs](http://www.popgen.dk/software/index.php/NgsAdmix))
    - `ibs:` Settings for identity by state calculation with ANGSD
      - `-doIBS:` Whether to use a random (1) or consensus (2) base in IBS
        distance calculation ([docs](http://www.popgen.dk/angsd/index.php/PCA_MDS))
