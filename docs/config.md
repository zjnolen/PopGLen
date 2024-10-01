# Configuring the workflow

Running the workflow requires configuring three files: `config.yaml`,
`samples.tsv`, and `units.tsv`. `config.yaml` is used to configure the
analyses, `samples.tsv` categorizes your samples into groups, and `units.tsv`
connects sample names to their input data files. The workflow will use
`config/config.yaml` automatically, but it can be good to name it something
informative and point to it when running snakemake with `--configfile <path>`.

### `samples.tsv`

This file contains your sample list, and has four **tab separated** columns:

```
sample	population	time	depth
hist1	Hjelmseryd	historical	low
hist2	Hjelmseryd	historical	low
hist3	Hjelmseryd	historical	low
mod1	Gotafors	modern	high
mod2	Gotafors	modern	high
mod3	Gotafors	modern	high
```

- `sample` contains the ID of a sample. It is best if it only contains
  alphanumeric characters.

- `population` contains the population the sample comes from and will be used
  to group samples for population-level analyses.

- `time` sets whether a sample should be treated as fresh DNA or historical DNA
  in the sequence processing workflow. Doesn't change anything if you're
  starting with bam files.

- `depth` puts the sample in a sequencing depth category. Used for filtering -
  if enabled in the configuration, extreme depth filters will be performed for
  depth categories individually.

### `units.tsv`

This file connects your samples to input files and has a potential for eight
**tab separated** columns:

```
sample	unit	lib	platform	fq1	fq2	bam	sra
hist1	BHVN22DSX2.2	hist1	ILLUMINA	data/fastq/hist1.r1.fastq.gz	data/fastq/hist1.r2.fastq.gz
hist1	BHVN22DSX2.3	hist1	ILLUMINA	data/fastq/hist1.unit2.r1.fastq.gz	data/fastq/hist1.unit2.r2.fastq.gz
hist2	BHVN22DSX2.2	hist2	ILLUMINA	data/fastq/hist2.r1.fastq.gz	data/fastq/hist2.r2.fastq.gz
hist3	BHVN22DSX2.2	hist2	ILLUMINA	data/fastq/hist3.r1.fastq.gz	data/fastq/hist3.r2.fastq.gz
mod1	AHW5NGDSX2.3	mod1	ILLUMINA	data/fastq/mod1.r1.fastq.gz	data/fastq/mod1.r2.fastq.gz
mod2	AHW5NGDSX2.3	mod2	ILLUMINA			data/bam/mod2.bam
mod3	AHW5NGDSX2.3	mod3	ILLUMINA	data/fastq/mod3.r1.fastq.gz	data/fastq/mod3.r2.fastq.gz
SAMN13218652	SRR10398077	SAMN13218652	ILLUMINA				SRR10398077
```

- `sample` contains the ID of a sample. Must be same as in `samples.tsv` and
  may be listed multiple times when inputting multiple sequencing
  runs/libraries.
- `unit` contains the sequencing unit, i.e. the sequencing lane barcode and
  lane number. This is used in the PU and (part of) the ID read groups. If you
  don't have multiple sequencing lanes per samples, this won't impact anything.
  Doesn't do anything when using bam input.
- `lib` contains the name of the library identifier for the entry. Fills in
  the LB and (part of) the ID read groups and is used for PCR duplicate removal.
  Best practice would be to have the combination of `unit` and `lib` to be
  unique per line. An easy way to use this is to use the Illumina library
  identifier or another unique library identifier, or simply combine a generic
  name with the sample name (sample1A, sample1B, etc.). Doesn't do anything when
  using bam input.
- `platform` is used to fill the PL read group. Commonly is just 'ILLUMINA'.
  Doesn't do anything when using bam input.
- `fq1` and `fq2` provides the absolute or relative to the working directory
  paths to the raw sequencing files corresponding to the metadata in the
  previous columns.
- `bam` provides the absolute or relative to the working directory path of
  pre-processed bam files. Only one bam files should be provided per sample in
  the units file.
- `sra` provides the NCBI SRA accession number for a set of paired end fastq
  files that will be downloaded to be processed. If a sample has multiple runs
  you would like to include, each run should be its own line in the units sheet,
  just as separate sequencing runs would be.

!!! note "Mixing samples with different starting points"
It is possible to have different samples start from different inputs (i.e.
some from bam, others from fastq, others from SRA). It is best to provide
only `fq1`+`fq2`, `bam`, or `sra` for a single sample to be clear where that
sample starts. If multiple starts are provided for the same sample, the bam file
will override fastq or SRA entries, and the fastq will override SRA
entries. Note that this means it is not currently possible to have multiple
starting points for _the same_ sample (i.e. FASTQ reads that would be
processed then merged into an existing BAM).

## Configuration file

`config.yaml` contains the configuration for the workflow, this is where you
will put what analyses, filters, and options you want. Below I describe the
configuration options. The [`config.yaml`](../config/config.yaml) in this
repository serves as a template and includes some 'default' parameters that may
be good starting points for many users. If `--configfile` is not specified in
the snakemake command, the workflow will default to `config/config.yaml`.

### Configuration options

#### Dataset Configuration

Required configuration of the 'dataset'.

- `samples:` An absolute or relative path from the working directory to the
  `samples.tsv` file.
- `units:` An absolute or relative path from the working directory to the
  `units.tsv` file.
- `dataset:` A name for this dataset run - essentially, an identifier for a
  batch of samples to be analysed together with the same configuration.

It is best to name your dataset something descriptive, but concise. This is
because this name will be used in organizing the results. Outputs of analyses
will be placed in the folder `results/{dataset}`, and files will be prefaced
with the dataset. This allows for multiple datasets to be run in the same
working directory, even in parallel (if they aren't trying to make the same
files), which is useful for multi-species projects or for testing out different
filters. You can simply have a config for each dataset and choose which one to
run with `--configfile`. A similar approach can be used to trying out different
analysis parameters.

??? note "Example use of multiple datasets
    Say you want to run PopGLen on two sets of samples, but use the same
    reference. You can have two sample lists: `dataset1_samples.tsv` and
    `dataset2_samples.tsv`, and two config files: `dataset1_config.tsv` and
    `dataset2_config.yaml`. In the configs, the `samples:` option should point
    to the corresponding sample list. The workflow for dataset1 can be run, if
    you pass `--configfile config/dataset1_config.yaml` to Snakemake, then the
    same can be done for dataset2. However, when dataset2 is run, it will use
    any outputs from dataset1 it can, such as reference indices, reference
    filters, etc. Additionally, if the two datasets share samples, those samples
    will not have to be remapped for dataset2, they'll be taken from the
    dataset1 run. The actual analyses are partitioned by dataset, into the
    folders `results/dataset1` and `results/dataset2`.

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

- `ancestral:` A path to a fasta file containing the ancestral states in your
  reference genome. This is optional, and is used to polarize allele
  frequencies in SAF files to ancestral/derived. If you leave this empty,
  the reference genome itself will be used as ancestral, and you should be
  sure the [`params`] [`realSFS`] [`fold`] is set to `1`. If you put a fasta
  here, you can set that to `0`.

Reference genomes should be uncompressed, and contig names should be clear and
concise. Currently, there are some issues parsing contig names with
underscores, so please change these in your reference before running the
pipeline. Alphanumeric characters, as well as `.` in contig names have been
tested to work so far, other symbols have not been tested.

Potentially the ability to use bgzipped genomes will be added, I just need to
check that it works with all underlying tools. Currently, it will for sure not
work, as calculating chunks is hard-coded to work on an uncompressed genome.

#### Sample Set Configuration

- `exclude_ind:` Sample name(s) that will be excluded from the workflow. Should
  be a list in []. Putting a `#` in front of the sample in the sample list also
  works. Mainly used to drop samples with poor quality after initial processing.
- `excl_pca-admix:` Sample name(s) that will be excluded _only_ from PCA and
  Admixture analyses. Useful for close relatives that violate the assumptions
  of these analyses, but that you want in others. Should be a list in []. If you
  want relatives out of all downstream analyses, not just PCA/Admix, put them in
  `exclude_ind` instead. Note this will trigger a re-run for relatedness
  analyses, but you can just disable them now as they've already been run.

#### Analysis Selection

Here, you will define which analyses you will perform. It is useful to start
with only a few, and add more in subsequent workflow runs, just to ensure you
catch errors before you use compute time running all analyses. Most are set
with (`true`/`false`) or a value, described below. Modifications to the
settings for each analysis are set in the next section.

- `populations:` A list of populations found in your sample list to limit
  population analyses to. Might be useful if you want to perform individual
  analyses on some samples but not include them in any population level
  analyses. Leave blank (`[]`) if you want population level analyses on all the
  populations defined in your `samples.tsv` file.

- `analyses:`
  - `mapping:`
    - `historical_only_collapsed:` Historical samples are expected to have
      fragmented DNA. For this reason, overlapping (i.e. shorter, usually
      <270bp) read pairs are collapsed in this workflow for historical samples.
      Setting this option to `true` will only map only these collapsed reads,
      and is recommended to target primarily endogenous content. However, in
      the event you want to map both the collapsed and uncollapsed reads, you
      can set this to `false`. (`true`/`false`)
    - `historical_collapsed_aligner:` Aligner used to map collapsed historical
      sample reads. `aln` is the recommended for this, but this is here in case
      you would like to select `mem` for this. Uncollapsed historical reads
      will be mapped with `mem` if `historical_only_collapsed` is set to
      `false`, regardless of what is put here. (`aln`/`mem`)
  - `pileup-mappability:` Filter out sites with low 'pileup mappability', which
    describes how uniquely fragments of a given size can map to the reference
    (`true`/`false`)
  - `repeatmasker:` (NOTE: Only one of the four options should be filled/true)
    - `bed:` Supply a path to a bed file that contains regions with repeats.
      This is for those who want to filter out repetitive content, but don't
      need to run Repeatmodeler or masker in the workflow because it has
      already been done for the genome you're using. Be sure the contig names
      in the bed file match those in the reference supplied. GFF or other
      filetypes that work with `bedtools subtract` may also work, but haven't
      been tested.
    - `local_lib:` Filter repeats by masking with an already made library you
      have locally (such as ones downloaded for Darwin Tree of Life genomes).
      Should be file path, not a URL.
    - `dfam_lib:` Filter repeats using a library available from dfam. Should be
      a taxonomic group name.
    - `build_lib:` Use RepeatModeler to build a library of repeats from the
      reference itself, then filter them from analysis (`true`/`false`).
  - `extreme_depth:` Filter out sites with extremely high or low global
    sequencing depth. Set the parameters for this filtering in the `params`
    section of the yaml. (`true`/`false`)
  - `dataset_missing_data:` A floating point value between 0 and 1. Sites with
    data for fewer than this proportion of individuals across the whole dataset
    will be filtered out in all analyses using the filtered sites file. (This is
    only needed if you need to ensure all your analyses are using exactly the
    same sites, which I find may result in coverage biases in results,
    especially heterozygosity. Unless you explicitly need to ensure all groups
    and analyses use the same sites, I would leave this blank, instead using
    the \[`params`]\[`angsd`]\[`minind_dataset`] to set a minimum individual
    threshold for dataset level analyses, allowing analyses to maximize sites
    per group/sample. This is how most papers do it.)
  - `population_missing_data:` A floating point value between 0 and 1. Sites
    with data for fewer than this proportion of individuals in any population
    will be filtered out in all populations using the filtered sites file.
    (This is only needed if you need to ensure all your populations are using
    exactly the same sites, which I find may result in coverage biases in
    results, especially heterozygosity. Unless you explicitly need to ensure all
    groups and analyses use the same sites, I would leave this blank, instead
    using the \[`params`]\[`angsd`]\[`minind_pop`] to set a minimum individual
    threshold for each analyses, allowing analyses to maximize sites per
    group/sample. This is how most papers do it.)
  - `qualimap:` Perform Qualimap bamqc on bam files for general quality stats
    (`true`/`false`)
  - `ibs_ref_bias:` Enable reference bias calculation. For each sample, one read
    is randomly sampled at each position and compared to the reference base.
    These are summarized as the proportion of the genome that is identical by
    state to the reference for each sample to quantify reference bias. This is
    done for all filter sets as well as for all sites without site filtering.
    If transition removal or other arguments are passed to ANGSD, they are
    included here. (`true`/`false`)
  - `damageprofiler:` Estimate post-mortem DNA damage on historical samples
    with Damageprofiler (`true`/`false`) NOTE: This just adds the addition of
    Damageprofiler to the already default output of MapDamage.
  - `mapdamage_rescale:` Rescale base quality scores using MapDamage2 to help
    account for post-mortem damage in analyses (if you only want to assess
    damage, use damageprofiler instead, they return the same results)
    (`true`/`false`) [docs](https://ginolhac.github.io/mapDamage/)
  - `estimate_ld:` Estimate pairwise linkage disquilibrium between sites with
    ngsLD for each popualation and the whole dataset. Note, only set this if
    you want to generate the LD estimates for use in downstream analyses
    outside this workflow. Other analyses within this workflow that require LD
    estimates (LD decay/pruning) will function properly regardless of the
    setting here. (`true`/`false`)
  - `ld_decay:` Use ngsLD to plot LD decay curves for each population and for
    the dataset as a whole (`true`/`false`)
  - `pca_pcangsd:` Perform Principal Component Analysis with PCAngsd. Currently
    requires at least 4 samples to finish, as it will by default try to plot
    PCs1-4. (`true`/`false`)
  - `admix_ngsadmix:` Perform admixture analysis with NGSadmix (`true`/`false`)
  - `relatedness:` Relatedness is estimated using two methods: IBSrelate (Waples
    et al. 2019, _Mol. Ecol._) and NgsRelate v2 (Hanghøj et al. 2019;
    _GigaScience_). IBSrelate does not require allele frequencies, which is
    useful if you do not have sufficient sample size to confidently estimate
    allele frequencies for your populations. In this pipeline, it is can be run
    three ways: using the (1) IBS and (2) SFS based methods described in the
    Waples paper using ANGSD or (3) using the SFS based method's implementation
    in NgsRelate v2 (which still does not require allele frequencies). NgsRelate
    v2 also offers an allele frequency based method, which enables co-inference
    of inbreeding and relatedness coefficients. If using this method, PopGLen
    will calculate the allele frequencies for your populations and input them
    into NgsRelate. These different methods have trade-offs in memory usage and
    run time. Generally, I recommend starting with the NgsRelate, using
    IBSrelate only (`ngsrelate_ibsrelate-only`), using the other approaches as
    you need them.
    - `ibsrelate_ibs:` Estimate pairwise relatedness with the IBS based method
      from Waples et al. 2019, _Mol. Ecol._. This can use a lot of memory, as
      it has genotype likelihoods for all sites from all samples loaded into
      memory, so it is done per 'chunk', which still takes a lot of time and
      memory. NOTE: For those removing transitions, this method does not include
      transition removal. All other relatedness methods here do.
      (`true`/`false`)
    - `ibsrelate_sfs:` Estimate pairwise relatedness with the SFS based method
      from Waples et al. 2019, _Mol. Ecol._. Enabling this can greatly increase
      the time needed to build the workflow DAG if you have many samples. As a
      form of this method is implemented in NGSrelate, it may be more
      efficient to only enable that. (`true`/`false`)
    - `ngsrelate_ibsrelate-only:` Performs the IBSrelate SFS method, but on SNPs
      using NgsRelate. Does not need to estimate allele frequencies.
      (`true`/`false`)
    - `ngsrelate_freqbased:` Performs the allele frequency based co-inference of
      relatedness and inbreeding that NgsRelate is primarly intended for. Will
      estimate allele frequencies per population and use them in the analysis.
      Also runs the IBSrelate SFS method in NgsRelate. (`true`/`false`)
  - `1dsfs:` Generates a one dimensional site frequency spectrum for all
    populations in the sample list. Automatically enabled if `thetas_angsd` is
    enabled. (`true`/`false`)
  - `1dsfs_boot:` Generates N bootstrap replicates of the 1D site frequency
    spectrum for each population. N is determined from the `sfsboot` setting
    below (`true`/`false`)
  - `2dsfs:` Generates a two dimensional site frequency spectrum for all unique
    populations pairings in the sample list. Automatically enabled if
    `fst_angsd` is enabled. (`true`/`false`)
  - `2dsfs_boot:` Generates N bootstrap replicates of the 2D site frequency
    spectrum for each population pair. N is determined from the `sfsboot`
    setting below (`true`/`false`)
  - `thetas_angsd:` Estimate pi, theta, and Tajima's D for each population in
    windows across the genome using ANGSD (`true`/`false`)
  - `heterozygosity_angsd:` Estimate individual genome-wide heterozygosity
    using ANGSD. Calculates confidence intervals from bootstraps.
    (`true`/`false`)
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
  - `pop_allele_freqs:` Estimates population-specific minor allele frequencies
    for each population in the dataset using ANGSD. Two outputs are generated
    per population: (1) population-specific minor allele frequencies, where only
    sites variable in the population are included and the minor allele is set
    to the minor of the population, and (2) dataset-wide minor allele
    frequencies, where the minor allele is set to the minor of the entire
    dataset and includes sites that are fixed within the population if they are
    variable in the dataset.

#### Subsampling Section

As this workflow is aimed at low coverage samples, its likely there might be
considerable variance in sample depth. For this reason, it may be good to
subsample all your samples to a similar depth to examine if variation in depth
is influencing results. To do this, set an integer value here to subsample all
your samples down to and run specific analyses. This subsampling can be done
in reference to the unfiltered sequencing depth, the mapping and base quality
filtered sequencing depth, or the filtered sites sequencing depth. The latter
is recommended, as this will ensure that sequencing depth is made uniform at
the analysis stage, as it is these filtered sites that analyses are performed
for.

- `subsample_dp:` A list of mean depths to subsample your reads to. This will be
  done per sample, and subsample from all the reads. Leaving list empty disables
  subsampling, list can contain any number of depths to subsample to. If a
  sample already has the same, or lower, depth than this number, it will just be
  used as is in the analysis. (List of INT)
- `subsample_by:` This determines how the 'full' sequencing depth of a sample
  is calculated to determine the amount of subsampling needed to reach the
  target depth. This should be one of three options: (1) `"unfilt"` will treat
  the sequencing depth of all (unfiltered) reads and sites as the 'full' depth;
  (2) `"mapqbaseq"` will filter out reads that don't pass the configured
  mapping or base quality, then calculate depth across all sites as the 'full'
  depth, (3) `"sitefilt"` will filter reads justa as `"mapqbaseq"` does, but
  will only calculate the 'full' depth from sites passing the sites filter. As
  the main goal of subsampling is to make depth uniform for analyses, this
  latter option is preferred, as it will most accurately bring the depth of the
  samples to the target depth for analyses.
  (`"unfilt"`/`"mapqbaseq"`/`"sitefilt"`)
- `redo_depth_filts`: If `subsample_by` is set to `"unfilt"` or `"mapqbaseq"`,
  then it would be possible to recaculate extreme depth filters for the
  subsampled dataset. Enable this to do so, otherwise, the depth filters from
  the full depth bams will be used. If `subsample_by` is set to `"sitefilt"`
  this will have no effect, as the subsampling is already in reference to a set
  site list. (`true`/`false`)
- `drop_samples`: When performing depth subsampling, you may want to leave some
  samples out that you kept in your 'full' dataset. These can be listed here and
  they will be removed from ALL depth subsampled analyses. A use case for this
  might be if you have a couple samples that are below your targeted subsample
  depth, and you don't want to include them. Note that if you configure multiple
  `subsample_dp`, these samples will be dropped from all of them. If you need to
  perform mutliple depth subsamplings with different subsets of samples, its
  best to run each depth individually. Alternatively, a config file can be made
  for each subsampled depth, however you may run into issues of file locking
  blocking both from running at the same time. (list of strings: `[]`)
- `subsample_analyses:` Individually enable analyses to be performed with the
  subsampled data. These are the same as the ones above in the analyses
  section. Enabling here will only run the analysis for the subsampled data,
  if you want to run it for the full data as well, you need to enable it in the
  analyses section as well. (`true`/`false`)

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
set of neutral sites you supply.

To skip running any analyses for `allsites-filts` and only perform them for the
BED files you supply, you can set `only_filter_beds: true` in the config file.
This may also be useful in the event you have a set of already filtered sites,
and want to run the workflow on those, ignoring any of the built in filter
options by setting them to `false`.

#### Software Configuration

This section contains the specific settings for each software, allowing users to
customize the settings used. The default configuration file contains settings
that are commonly used, and should be applicable to most datasets sequenced on
patterened flow cells, but please check that they make sense for your analysis.
If you are missing a configurable setting you need, open up an issue or a pull
request and I'll gladly put it in if possible.

!!! note "Note to historical sample users wanting to remove transitions"
    While most the defaults below are good for most datasets, including ones
    with historical samples and using MapDamage rescaling, transition removal is
    turned off by default. To enable transition removal to account for
    post-mortem DNA damage, enable the option `rmtrans` in the `angsd` section
    below. This will fill in the appropriate flag `-rmTrans` or `-noTrans`
    depending on the analysis, and remove transitions from all analyses. Only
    the IBS based IBSrelate method currently does not support transition
    removal.

- `mapQ:` Phred-scaled mapping quality filter. Reads below this threshold will
  be filtered out. (integer)
- `baseQ:` Phred-scaled base quality filter. Reads below this threshold will be
  filtered out. (integer)

- `params:`
  - `clipoverlap:`
    - `clip_user_provided_bams:` Determines whether overlapping read pairs will
      be clipped in BAM files supplied by users. This is useful as many variant
      callers will account for overlapping reads in their processing, but ANGSD
      will double count overlapping reads. If BAMs were prepped without this in
      mind, it can be good to apply before running through ANGSD. However, it
      essentially creates a BAM file of nearly equal size for every sample, so
      it may be nice to turn off if you don't care for this correction or have
      already applied it on the BAMs you supply. (`true`/`false`)
  - `genmap:` Parameters for pileup mappability analysis, see
    [GenMap's documentation](https://github.com/cpockrandt/genmap/) for more
    details.
    - `K:`
    - `E:`
    - `map_thresh:` A threshold mappability score. Each site gets an average
      mappability score taken by averaging the mappability of all K-mers that
      would overlap it. A score of 1 means all K-mers are uniquely mappable,
      allowing for `e` mismatches. This is doen via a custom script, and may
      eventually be replaced by the SNPable method, which is more common.
      (integer/float, 0-1)
  - `extreme_depth_filt:` Parameters for excluding sites based on extreme high
    and/or low global depth. The final sites list will contain only sites that
    pass the filters for all categories requested (i.e the whole dataset
    and/or the depth categories set in samples.tsv).
    - `method:` Whether you will generate extreme thresholds as a multiple of
      the median global depth (`"median"`) or as percentiles of the
      global depth distribution (`"percentile"`)
    - `bounds:` The bounds of the depth cutoff, defined as a numeric list. For
      the median method, the values will be multiplied by the median of the
      distribution to set the thresholds (i.e. `[0.5,1.5]` would generate
      a lower threshold at 0.5\*median and an upper at 1.5\*median). For the
      percentile method, these define the lower and upper percentiles to filter
      out (i.e [0.01,0.99] would remove the lower and upper 1% of the depth
      distributions). (`[ FLOAT, FLOAT]`)
    - `filt-on-dataset:` Whether to perform this filter on the dataset as a
      whole (may want to set to false if your dataset global depth distribution
      is multi-modal). (`true`/`false`)
    - `filt-on-depth-classes:` Whether to perform this filter on the depth
      classes defined in the samples.tsv file. This will generate a global
      depth distribution for samples in the same category, and perform the
      filtering on these distributions independently. Then, the sites that pass
      for all the classes will be included. (`true`/`false`)
  - `fastp:`
    - `extra:` Additional options to pass to fastp trimming. (string)
    - `min_overlap_hist:` Minimum overlap to collapse historical reads. Default
      in fastp is 30. This effectively overrides the `--length_required` option
      if it is larger than that. (INT)
  - `bwa_aln:`
    - `extra:` Additional options to pass to bwa aln for mapping of historical
      sample reads. (string)
  - `samtools:`
    - `subsampling_seed:` Seed to use when subsampling bams to lower depth.
      `"$RANDOM"` can be used to set a random seed, or any integer can be used
      to set a consistent seed. (string or int)
  - `picard:`
    - `MarkDuplicates:` Additional options to pass to Picard MarkDuplicates.
      `--REMOVE_DUPLICATES true` is recommended. (string)
  - `angsd:` General options in ANGSD, relevant doc pages are linked
    - `gl_model:` Genotype likelihood model to use in calculation
      (`-GL` option in ANGSD, [docs](http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods))
    - `maxdepth:` When calculating individual depth, sites with depth higher
      than this will be binned to this value. Should be fine for most to leave
      at `1000`. (integer, [docs](http://www.popgen.dk/angsd/index.php/Depth))
    - `mindepthind:` Individuals with sequencing depth below this value at a
      position will be treated as having no data at that position by ANGSD.
      ANGSD defaults to 1 for this. Note that this can be separately set for
      individual heterozygosity estimates with `mindepthind_heterozygosity`
      below. (integer, `-setMinDepthInd` option in ANGSD) (INT)
    - `minind_dataset:` Used to fill the `-minInd` option for any dataset wide
      ANGSD outputs (like Beagles for PCA/Admix). Should be a floating point
      value between 0 and 1 describing what proportion of the dataset must have
      data at a site to include it in the output. (FLOAT)
    - `minind_pop:` Used to fill the `-minInd` option for any population level
      ANGSD outputs (like SAFs or Beagles for ngsF-HMM). Should be a floating
      point value between 0 and 1 describing what proportion of the population
      must have data at a site to include it in the output. (FLOAT)
    - `rmtrans:` Removes transitions using ANGSD, effectively removing them
      from downstream analyses. This is useful for removing DNA damage from
      analyses, and will automatically set the appropriate ANGSD flags (i.e.
      using `-noTrans 1` for SAF files and `-rmTrans 1` for Beagle files.)
      (`true`/`false`)
    - `extra:` Additional options to pass to ANGSD during genotype likelihood
      calculation at all times. This is primarily useful for adding BAM input
      filters. Note that `--remove_bads` and `-only_proper_pairs` are enabled
      by default, so they only need to be included if you want to turn them
      off or explicitly ensure they are enabled. I've also found that for some
      datasets, `-C 50` and `-baq 1` can create a strong relationship between
      sample depth and detected diversity, effectively removing the benefits of
      ANGSD for low/variable depth data. I recommend that these aren't included
      unless you know you need them. Since the workflow uses bwa to map,
      `-uniqueOnly 1` doesn't do anything if your minimum mapping quality is
      \> 0. Mapping and base quality thresholds are also not needed, it will
      use the ones defined above automatically. If you prefer to correct for
      historical damage by trimming the ends of reads, this is where you'd want
      to put `-trim INT`. (string)
      (string, [docs](http://www.popgen.dk/angsd/index.php/Input#BAM.2FCRAM))
    - `extra_saf:` Same as `extra`, but only used when making SAF files (used
      for SFS, thetas, Fst, IBSrelate, heterozygosity includes invariable
      sites). Doesn't require options already in `extra` or defined via other
      params in the YAML (such as `notrans`, `minind`, `GL`, etc.) (string)
    - `extra_beagle:` Same as `extra`, but only used when making Beagle and Maf
      files (used for PCA, Admix, ngsF-HMM, doIBS, ngsrelate, includes only
      variable sites). Doesn't require options already in `extra` or defined via
      other params in the YAML (such as `rmtrans`, `minind`, `GL`, etc.)
      (string)
    - `snp_pval:` The p-value to use for calling SNPs
      (float, [docs](http://www.popgen.dk/angsd/index.php/SNP_calling)) (float
      or string)
    - `domajorminor:` Method for inferring the major and minor alleles. Set to
      1 to infer from the genotype likelihoods, see
      [documentation](https://www.popgen.dk/angsd/index.php/Major_Minor)
      for other options. `1`, `2`, and `4` can be set without any additional
      configuration. `5` must also have an ancestral reference provided in the
      config, otherwise it will be the same as `4`. `3` is currently not
      possible, and is used for generating the dataset minor allele MAFs. If you
      have a use for `3` that PopGLen is suited for, please open an issue and I
      will look into it. (int)
    - `domaf:` Method for inferring minor allele frequencies. Set to `1` to
      infer from genotype likelihoods using a known major and minor from the
      `domajorminor` setting above. Set to `2` to assume the minor is unknown.
      See [docs](http://www.popgen.dk/angsd/index.php/Allele_Frequencies) for
      more information. `3` is possible, and will estimate the frequencies both
      assuming a known and unknown minor. If you choose this option, you'll get
      both in the MAF outputs, but only the known will be passed to NgsRelate if
      you also use that. Other values are currently unsupported in PopGLen.
      (int)
    - `min_maf:` The minimum minor allele frequency required to call a SNP.
      This is set when generating the beagle file, so will filter SNPs for
      PCAngsd, NGSadmix, ngsF-HMM, and NGSrelate. If you would like each tool
      to handle filtering for maf on its own you can set this to `-1`
      (disabled). (float, [docs](http://www.popgen.dk/angsd/index.php/Allele_Frequencies))
    - `mindepthind_heterozygosity:` When estimating individual heterozygosity,
      sites with sequencing depth lower than this value will be dropped.
      (integer, `-setMinDepthInd` option in ANGSD) (int)
  - `ngsld:` Settings for ngsLD ([docs](https://github.com/fgvieira/ngsLD))
    - `max_kb_dist_est-ld:` For the LD estimates generated when setting
      `estimate_ld: true` above, set the maximum distance between sites in kb
      that LD will be estimated for (`--max_kb_dist` in ngsLD, integer)
    - `rnd_sample_est-ld:` For the LD estimates generated when setting
      `estimate_ld: true` above, randomly sample this proportion of pairwise
      linkage estimates rather than estimating all (`--rnd_sample` in ngsLD,
      float)
    - `max_kb_dist_decay:` The same as `max_kb_dist_est-ld:`, but used when
      estimating LD decay when setting `ld_decay: true` above (integer)
    - `rnd_sample_decay:` The same as `rnd_sample_est-ld:`, but used when
      estimating LD decay when setting `ld_decay: true` above (float)
    - `fit_LDdecay_extra:` Additional plotting arguments to pass to
      `fit_LDdecay.R` when estimating LD decay (string)
    - `fit_LDdecay_n_correction:` When estimating LD decay, should the sample
      size corrected r^2 model be used? (`true`/`false`, `true` is the
      equivalent of passing a sample size to `fit_LDdecay.R` in ngsLD using
      `--n_ind`)
    - `max_kb_dist_pruning_dataset:` The same as `max_kb_dist_est-ld:`, but
      used when linkage pruning SNPs as inputs for PCAngsd, NGSadmix, and
      NGSrelate analyses. Pruning is performed on the whole dataset. Any
      positions above this distance will be assumed to be in linkage
      equilibrium during the pruning process. (integer)
    - `pruning_min-weight_dataset:` The minimum r^2 to assume two positions are
      in linkage disequilibrium when pruning for PCAngsd, NGSadmix, and
      NGSrelate analyses. (float)
  - `ngsrelate:` Settings for NGSrelate
    - `ibsrelate-only-extra:` Any extra command line parameters to be passed to
      NgsRelate when performing an IBSrelate only (no allele frequencies) run.
      (string)
    - `freqbased-extra:` Any extra command line parameters to be passed to
      NgsRelate when performing a standard, allele frequency based, run.
      (string)
  - `ngsf-hmm:` Settings for ngsF-HMM
    - `estimate_in_pops:` Set to `true` to run ngsF-HMM separately for each
      population in your dataset. Set to `false` to run for whole dataset at
      once. ngsF-HMM assumes Hardy-Weinberg Equilibrium (aside from inbreeding)
      in the input data, so select the option that most reflects this in your
      data. (`true`/`false`)
    - `prune:` Whether or not to prune SNPs for LD before running the analysis.
      ngsF-HMM assumes independent sites, so it is preferred to set this to
      `true` to satisfy that expectation. (`true`/`false`)
    - `max_kb_dist_pruning_pop:` The maximum distance between sites in kb
      that will be treated as in LD when pruning for the ngsF-HMM input. (INT)
    - `pruning_min-weight_pop:` The minimum r^2 to assume two positions are in
      linkage disequilibrium when pruning for the ngsF-HMM input. Note, that
      this likely will be substantially higher for individual populations than
      for the whole dataset, as background LD should be higher when no
      substructure is present. (float)
    - `min_roh_length:` Minimum ROH size in base pairs to include in inbreeding
      coefficient calculation. Set if short ROH might be considered low
      confidence for your data. (integer)
    - `roh_bins:` A list of integers that describe the size classes in base
      pairs you would like to partition the inbreeding coefficient by. This can
      help visualize how much of the coefficient comes from ROH of certain size
      classes (and thus, ages). List should be in ascending order and the first
      entry should be greater than `min_roh_length`. The first bin will group
      ROH between `min_roh_length` and the first entry, subsequent bins will
      group ROH with sizes between adjacent entries in the list, and the final
      bin will group all ROH larger than the final entry in the list. (list)
  - `realSFS:` Settings for realSFS
    - `fold:` Whether or not to fold the produced SFS. Set to 1 if you have not
      provided an ancestral-state reference (0 or 1,
      [docs](http://www.popgen.dk/angsd/index.php/SFS_Estimation))
    - `sfsboot:` Determines number of bootstrap replicates to use when
      requesting bootstrapped SFS. Is used for both 1dsfs and 2dsfs (this is
      very easy to separate, open an issue if desired). Automatically used
      for heterozygosity analysis to calculate confidence intervals. (integer)
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
    - `minsites:` Minimum sites to include window in report plot. This does not
      remove them from the actual output, just the report plot.
  - `ngsadmix:` Settings for admixture analysis with NGSadmix. This analysis is
    performed for a set of K groupings, and each K has several replicates
    performed. Replicates will continue until a set of N highest likelihood
    replicates converge, or the number of replicates reaches an upper threshold
    set here. Defaults for `reps`, `minreps`, `thresh`, and `conv` can be left
    as default for most.
    - `kvalues:` A list of values of K to fit the data to (list of integers)
    - `reps:` The maximum number of replicates to perform per K. Default is 100.
      (integer)
    - `minreps:` The minimum number of replicates to perform, even if
      replicates have converged. Default is 20. (integer)
    - `thresh:` The convergence threshold - the top replicates must all be
      within this value of log-likelihood units to consider the run converged.
      Default is 2. (integer)
    - `conv:` The number of top replicates to include in convergence
      assessment. Default is 3. (integer)
    - `extra:` Additional arguments to pass to NGSadmix (for instance,
      increasing `-maxiter`). (string,
      [docs](http://www.popgen.dk/software/index.php/NgsAdmix))
  - `ibs:` Settings for identity by state calculation with ANGSD
    - `-doIBS:` Whether to use a random (1) or consensus (2) base in IBS
      distance calculation
      ([docs](http://www.popgen.dk/angsd/index.php/PCA_MDS))
