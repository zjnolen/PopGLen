# Configuration cheat sheet

This contains the same information largely as the
[configuration docs](config.md), but in an interactive `config.yaml` file. This
can be helpful if you have your config in front of you and want to check what a
setting does. Just click the `+` next to each line to get info on that setting.

```yaml title="config/config.yaml" linenums="1"
#=====================Dataset Configuration============================#

samples: config/samples.tsv # (1)

units: config/units.tsv # (2)

dataset: "dataset-name" # (3)

#===================== Reference Configuration ========================#

chunk_size: 50000000 # (4)

reference:
  name: "dataset-ref" # (5)
  fasta: "resources/ref/reference-genome.fa" # (6)
  mito: ["mitochondrion"] # (7)
  sex-linked: ["Z"] # (8)
  exclude: [] # (9)
  min_size: 1000000 # (10)

ancestral: # (11)

#===================== Sample Set Configuration =======================#

exclude_ind: [] # (12)

excl_pca-admix: [] # (13)

#====================== Analysis Selection ============================#

populations: [] # (14)

analyses:
  # mapping options
  mapping:
    historical_only_collapsed: true # (15)
    historical_collapsed_aligner: "aln" # (16)
  # sites file filters
  pileup-mappability: # (17)
  repeatmasker: # (18)
    bed: # (19)
    local_lib: # (20)
    build_lib: true # (21)
  extreme_depth: true # (22)
  dataset_missing_data: # (23)
  population_missing_data: # (24)
  # quality control
  qualimap: true # (25)
  ibs_ref_bias: true # (26)
  damageprofiler: true # (27)
  mapdamage_rescale: # (28)
  # population genomic analyses
  estimate_ld: # (29)
  ld_decay: true # (30)
  pca_pcangsd: true # (31)
  admix_ngsadmix: true # (32)
  relatedness: # (33)
    ibsrelate_ibs: # (34)
    ibsrelate_sfs: # (35)
    ngsrelate_ibsrelate-only: true # (36)
    ngsrelate_freqbased: # (37)
  1dsfs: # (38)
  1dsfs_boot: # (39)
  2dsfs: # (40)
  2dsfs_boot: # (41)
  thetas_angsd: true # (42)
  heterozygosity_angsd: true # (43)
  fst_angsd: # (44)
    populations: true # (45)
    individuals: # (46)
  inbreeding_ngsf-hmm: true # (47)
  ibs_matrix: # (48)
  pop_allele_freqs: true # (49)

#==================== Downsampling Configuration ======================#

subsample_dp: [4] # (50)

subsample_by: "sitefilt" # (51) three options: unfilt, mapqbaseq, sitefilt; sitefilt recommended (subsamples bams to `subsample_dp` within filtered regions)
redo_depth_filts: false # (52) has no effect when `subsample_by` == "sitefilt"
drop_samples: [] # (53) These samples won't be included in the downsampled dataset

subsample_analyses: # (54)
  estimate_ld:
  ld_decay:
  pca_pcangsd:
  admix_ngsadmix:
  relatedness:
    ibsrelate_ibs:
    ibsrelate_sfs:
    ngsrelate_ibsrelate-only:
    ngsrelate_freqbased:
  1dsfs:
  1dsfs_boot:
  2dsfs:
  2dsfs_boot:
  thetas_angsd: true
  heterozygosity_angsd: true
  fst_angsd:
    populations: true
    individuals:
  inbreeding_ngsf-hmm: true
  ibs_matrix:
  pop_allele_freqs: true

#=========================== Filter Sets ==============================#

filter_beds: #(55)
  example:

only_filter_beds: false # (56)

#===================== Software Configuration =========================#

mapQ: 30 # (57)
baseQ: 20 # (58)

params:
  clipoverlap:
    clip_user_provided_bams: false # (59)
  fastp:
    extra: "-p -g" # (60) don't put --merge or --overlap_len_require here, they're implicit
    min_overlap_hist: 30 # (61)
  bwa_aln:
    extra: "-l 16500 -n 0.01 -o 2" # (62)
  picard:
    MarkDuplicates: "--REMOVE_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500" #(63)
  damageprofiler:
    profile_modern: false # (64)
  genmap: # (65)
    K: "30" # (66)
    E: "2" # (67)
    map_thresh: 1 # (68)
  extreme_depth_filt: # (69)
    method: "percentile" # (70)
    bounds: [ 0.025, 0.975 ] # (71)
    filt-on-dataset: true # (72)
    filt-on-depth-classes: true # (73)
  samtools:
    subsampling_seed: "$RANDOM" # (74)
  angsd: # (75)
    gl_model: 2 # (76) gl_model - 1=SAMtools, 2=GATK, 3=SOAPsnp, 4=SYK
    maxdepth: 100 # (77) for depth calc only
    mindepthind: 1 # (78)
    minind_dataset: 0.75 # (79) as a proportion of N
    minind_pop: 0.75 # (80) as a proportion of N
    rmtrans: true # (81)
    extra: "" # (82) goes to all ANGSD runs
    extra_saf: "" # (83) goes to all -doSaf runs
    extra_beagle: "" # (84) goes to all -doGlf 2 runs
    domajorminor: 1 # (85) Currently only compatible with values 1, 2, 4, 5
    domaf: 1 # (87) Currently only compatible with values 1 & 2
    snp_pval: "1e-6" # (86)
    min_maf: 0.05 # (88) Set -1 to disable
    mindepthind_heterozygosity: 1 # (89)
  ngsld: # (90)
    max_kb_dist_est-ld: 4000 # (91)
    rnd_sample_est-ld: 0.001 # (92)
    max_kb_dist_decay: 100 # (93)
    rnd_sample_decay: 0.01 # (94)
    fit_LDdecay_extra: "--fit_level 1 --plot_size 3,6 --plot_y_lim 1" # (95)
    fit_LDdecay_n_correction: true # (96)
    max_kb_dist_pruning_dataset: 25 # (97) used for PCA, Admix, not ngsF-HMM
    pruning_min-weight_dataset: 0.2 # (98) used for PCA, Admix, not ngsF-HMM
  ngsrelate: # (99)
    ibsrelate-only-extra: "" #(100)
    freqbased-extra: "" # (101)
  ngsf-hmm: # (102)
    estimate_in_pops: true # (103)
    prune: true # (104)
    max_kb_dist_pruning_pop: 25 # (105) used only for ngsf-hmm
    pruning_min-weight_pop: 0.4 # (106) used only for ngsf-hmm
    min_roh_length: 50000 # (107)
    roh_bins: [ 100000, 250000, 500000, 1000000, 2500000, 5000000 ] # (108) in ascending order, bins are between adjacent values and between min_roh_length and value 1 and the last value and infinity
  realsfs: # (109)
    fold: 1 # (110) Should only be set to 1, unless ancestral reference is given above
    sfsboot: 50 # (111)
  fst: # (112)
    whichFst: 1 # (113)
    win_size: 25000 # (114)
    win_step: 10000 # (115)
  thetas: # (116)
    win_size: 25000 # (117)
    win_step: 10000 # (118)
    minsites: 1000 # (119)
  ngsadmix: # (120)
    kvalues: [1,2,3,4,5] # (121)
    reps: 100 # (122)
    minreps: 20 # (123)
    thresh: 2 # (124)
    conv: 3 # (125)
    extra: "-maxiter 4000" # (126)
  ibs: # (127)
    doibs: 1 # (128)
```

1. An absolute or relative path from the working directory to the `samples.tsv`
   file.
2. An absolute or relative path from the working directory to the `units.tsv`
   file.
3. A name for this dataset run - an identifier for a batch of samples to be
   analysed together with the same configuration.
4. A size in bp (integer). Your reference will be analyzed in
   'chunks' of contigs of this size to parallelize processing. This size should
   be larger than the largest contig in your genome. A larger number means fewer
   jobs that run longer. A smaller number means more jobs that run shorter. The
   best fit will depend on the reference and the compute resources you have
   available. Leaving this blank will not divide the reference up into chunks.
5. A name for your reference genome, will go in the file names and be
   used to organize the reference filter outputs.
6. A path to the reference fasta file (currently only supports
   uncompressed fasta files).
7. Mitochondrial contig name(s), will be removed from analysis. Should
   be listed as a list of strings within brackets `[]`.
8. Sex-linked contig name(s), will be removed from analysis.
   Should be listed as a list of strings within brackets `[]`.
9. Additional contig name(s) to exclude from analysis. Should be
   listed as a list of strings within brackets `[]`.
10. A size in bp (integer). All contigs below this size will be excluded from
    analysis.
11. A path to a fasta file containing the ancestral states in your
    reference genome. This is optional, and is used to polarize allele
    frequencies in SAF files to ancestral/derived. If you leave this empty,
    the reference genome itself will be used as ancestral, and you should be
    sure the [`params`] [`realSFS`] [`fold`] is set to `1`. If you put a fasta
    here, you can set that to `0`.
12. Sample name(s) that will be excluded from the workflow. Should
    be a list in `[]`. As an alternative, putting a `#` in front of the sample
    in the sample list also will remove it. Mainly used to drop samples with
    poor quality after initial processing.
13. Sample name(s) that will be excluded _only_ from PCA and
    Admixture analyses. Useful for close relatives that violate the assumptions
    of these analyses, but that you want in others. Should be a list in `[]`. If
    you   want relatives out of all downstream analyses, not just PCA/Admix, put
    them in `exclude_ind` instead. Note this will trigger a re-run for
    relatedness analyses, but you can just disable them now as they've already
    been run.
14. A list of populations found in your sample list to limit
    population analyses to. Might be useful if you want to perform individual
    analyses on some samples, but not include them in any population level
    analyses. Leave blank (`[]`) if you want population level analyses on all
    the populations defined in your `samples.tsv` file.
15. Historical samples are expected to have
    fragmented DNA. For this reason, overlapping (i.e. shorter, usually
    <270bp) read pairs are collapsed in this workflow for historical samples.
    Setting this option to `true` will only map only these collapsed reads,
    and is recommended to target primarily endogenous content. However, in
    the event you want to map both the collapsed and uncollapsed reads, you
    can set this to `false`. (`true`/`false`)
16. Aligner used to map collapsed historical
    sample reads. `aln` is the recommended for this, but this is here in case
    you would like to select `mem` instead. Uncollapsed historical reads
    will be mapped with `mem` if `historical_only_collapsed` is set to
    `false`, regardless of what is put here. (`aln`/`mem`)
17. Filter out sites with low 'pileup mappability', which
    describes how uniquely fragments of a given size can map to the reference
    (`true`/`false`)
18. (NOTE: Only one of the three options should be filled/true)
19. Supply a path to a bed file that contains regions with repeats.
    This is for those who want to filter out repetitive content, but don't
    need to run Repeatmodeler or masker in the workflow because it has
    already been done for the genome you're using. Be sure the contig names
    in the bed file match those in the reference supplied. GFF or other
    filetypes that work with `bedtools subtract` may also work, but haven't
    been tested.
20. Filter repeats by supplying a repeat library you have locally
    (such as ones downloaded for Darwin Tree of Life genomes) and using it to
    identify repeats with RepeatMasker. Should be file path, not a URL.
21. Use RepeatModeler to build a library of repeats from the
    reference itself, then identify them with RepeatMasker and filter them out
    (`true`/`false`).
22. Filter out sites with extremely high or low global
    sequencing depth. Set the parameters for this filtering in the `params`
    section of the yaml. (`true`/`false`)
23. A floating point value between 0 and 1. Sites with
    data for fewer than this proportion of individuals across the whole dataset
    will be filtered out in all analyses using the filtered sites file. (This is
    only needed if you need to ensure all your analyses are using exactly the
    same sites, which I find may result in coverage biases in results,
    especially heterozygosity. Unless you explicitly need to ensure all groups
    and analyses use the same sites, I would leave this blank, instead using
    the \[`params`]\[`angsd`]\[`minind_dataset`] to set a minimum individual
    threshold for dataset level analyses, allowing analyses to maximize sites
    per group/sample. This is how most papers do it.)
24. A floating point value between 0 and 1. Sites
    with data for fewer than this proportion of individuals in any population
    will be filtered out in all populations using the filtered sites file.
    (This is only needed if you need to ensure all your populations are using
    exactly the same sites, which I find may result in coverage biases in
    results, especially heterozygosity. Unless you explicitly need to ensure all
    groups and analyses use the same sites, I would leave this blank, instead
    using the \[`params`]\[`angsd`]\[`minind_pop`] to set a minimum individual
    threshold for each analyses, allowing analyses to maximize sites per
    group/sample. This is how most papers do it.)
25. Perform Qualimap bamqc on bam files for general quality stats
    (`true`/`false`)
26. Enable reference bias calculation. For each sample, one read
    is randomly sampled at each position and compared to the reference base.
    These are summarized as the proportion of the genome that is identical by
    state to the reference for each sample to quantify reference bias. This is
    done for all filter sets as well as for all sites without site filtering.
    If transition removal or other arguments are passed to ANGSD, they are
    included here. (`true`/`false`)
27. Estimate post-mortem DNA damage on historical samples
    with Damageprofiler (`true`/`false`)
28. Rescale base quality scores using MapDamage2 to help
    account for post-mortem damage in analyses (if you only want to assess
    damage, use damageprofiler instead, they return the same results)
    (`true`/`false`) [[docs]](https://ginolhac.github.io/mapDamage/)
29. Estimate pairwise linkage disquilibrium between sites with
    ngsLD for each popualation and the whole dataset. Note, only set this if
    you want to generate the LD estimates for use in downstream analyses
    outside this workflow. Other analyses within this workflow that require LD
    estimates (LD decay/pruning) will function properly regardless of the
    setting here. (`true`/`false`) [[docs]](https://github.com/fgvieira/ngsLD)
30. Use ngsLD to plot LD decay curves for each population and for
    the dataset as a whole (`true`/`false`)
    [[docs]](https://github.com/fgvieira/ngsLD)
31. Perform Principal Component Analysis with PCAngsd. Currently
    requires at least 4 samples to finish, as it will by default try to plot
    PCs1-4. (`true`/`false`)
    [[docs]](https://popgen.dk/software/index.php/PCAngsd)
32. Perform admixture analysis with NGSadmix (`true`/`false`)
    [[docs]](https://www.popgen.dk/software/index.php/NgsAdmix)
33. Relatedness is estimated using two methods: IBSrelate and
    NgsRelate v2. IBSrelate does not require allele frequencies, which is
    useful if you do not have sufficient sample size to confidently estimate
    allele frequencies for your populations. In this pipeline, it is can be run
    three ways: using the (1) IBS and (2) SFS based methods described in Waples
    et al. 2019, _Mol. Ecol._ using ANGSD or (3) using the SFS based method's
    implementation in NgsRelate v2 (which still does not require allele
    frequencies). NgsRelate v2 also offers an allele frequency based method,
    which enables co-inference of inbreeding and relatedness coefficients. If
    using this method, PopGLen will calculate the allele frequencies for your
    populations and input them into NgsRelate. These different methods have
    trade-offs in memory usage and run time in addition to the methodological
    effects on the estimates themselves. Generally, I recommend starting with
    NgsRelate, using IBSrelate only (`ngsrelate_ibsrelate-only`), adding the
    other approaches as you need them.
34. Estimate pairwise relatedness with the IBS based method
    from. This can use a lot of memory, as it has genotype likelihoods for all
    sites from all samples loaded into memory, so it is done per 'chunk',
    which still takes a lot of time and memory. NOTE: For those removing
    transitions, this method does not include transition removal. All other
    relatedness methods here do. (`true`/`false`)
    [[docs]](https://popgen.dk/software/index.php/IBSrelate)
35. Estimate pairwise relatedness with the SFS based method.
    Enabling this can greatly increase the time needed to build the workflow
    DAG if you have many samples. As a form of this method is implemented in
    NGSrelate, it may be more efficient to only enable that. (`true`/`false`)
    [[docs]](https://popgen.dk/software/index.php/IBSrelate)
36. Performs the IBSrelate SFS method, but on SNPs
    only using NgsRelate. Does not need to estimate allele frequencies.
    (`true`/`false`) [[docs]](https://github.com/ANGSD/NgsRelate)
37. Performs the allele frequency based co-inference of
    relatedness and inbreeding that NgsRelate is primarly intended for. Will
    estimate allele frequencies per population and use them in the analysis.
    Also runs the IBSrelate SFS method in NgsRelate. (`true`/`false`)
    [[docs]](https://github.com/ANGSD/NgsRelate)
38. Generates a one dimensional site frequency spectrum for all
    populations in the sample list. Automatically enabled if `thetas_angsd` is
    enabled. (`true`/`false`)
    [[docs]](https://popgen.dk/angsd/index.php/SFS_Estimation)
39. Generates N bootstrap replicates of the 1D site frequency
    spectrum for each population. N is determined from the `sfsboot` setting
    below (`true`/`false`)
    [[docs]](https://popgen.dk/angsd/index.php/SFS_Estimation)
40. Generates a two dimensional site frequency spectrum for all unique
    populations pairings in the sample list. Automatically enabled if
    `fst_angsd` is enabled. (`true`/`false`)
    [[docs]](https://popgen.dk/angsd/index.php/SFS_Estimation)
41. Generates N bootstrap replicates of the 2D site frequency
    spectrum for each population pair. N is determined from the `sfsboot`
    setting below (`true`/`false`)
    [[docs]](https://popgen.dk/angsd/index.php/SFS_Estimation)
42. Estimate nucleotide diversity (pi), Watterson's theta, and
    Tajima's D for each population in windows across the genome using ANGSD. If
    polarized to an ancestral reference, additional measures of theta and
    neutrality statistics in the output will be relevant (see 'Unknown ancestral
    state (folded sfs)' in the ANGSD docs for more details on this)
    (`true`/`false`)
    [[docs]](https://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests)
43. Estimate individual genome-wide heterozygosity
    using ANGSD. Calculates confidence intervals from bootstraps.
    (`true`/`false`)
    [[docs]](https://www.popgen.dk/angsd/index.php/Heterozygosity)
44. Estimate pairwise _F_~ST~ using ANGSD. Set one or both of the
    below options. Estimates both globally and in windows across the genome.
45. Pairwise _F_~ST~ is calculated between all possible
    population pairs (`true`/`false`)
    [[docs]](https://popgen.dk/angsd/index.php/Fst)
46. Pairwise _F_~ST~ is calculated between all possible
    individual pairs. NOTE: This can be really intensive on the DAG building
    process, so I don't recommend enabling unless you're certain you want
    this (`true`/`false`) [[docs]](https://popgen.dk/angsd/index.php/Fst)
47. Estimates inbreeding coefficients and runs of
    homozygosity using ngsF-HMM. Inbreeding coefficients are output in both the
    model based form from ngsF-HMM, as well as converted into _F~RoH~_, which
    describes the proportion of the genome in runs of homozygosity over a
    certain length. (`true`/`false`)
    [[docs]](https://github.com/fgvieira/ngsF-HMM)
48. Estimate pairwise identity by state distance between all
    samples using ANGSD. (`true`/`false`)
    [[docs]](https://popgen.dk/angsd/index.php/PCA_MDS)
49. Estimates population-specific minor allele frequencies
    for each population in the dataset using ANGSD. Two outputs are generated
    per population: (1) population-specific minor allele frequencies, where only
    sites variable in the population are included and the minor allele is set
    to the minor of the population, and (2) dataset-wide minor allele
    frequencies, where the minor allele is set to the minor of the entire
    dataset and includes sites that are fixed within the population if they are
    variable in the dataset. Alternatively, major alleles can be polarized to
    the reference or ancestral state in both outputs if `domajorminor` is set to
    `4` or `5`. [[docs]](https://popgen.dk/angsd/index.php/Allele_Frequencies)
50. A list of mean depths to subsample your reads to. This will be
    done per sample, and subsample from all the reads. Leaving the list empty
    disables subsampling. The list can contain multiple depths to subsample to.
    If a sample already has the same, or lower, depth than this number, it will
    just be used as is in the analysis. (List `[]` of integers)
51. This determines how the 'full' sequencing depth of a sample
    is calculated to determine the amount of subsampling needed to reach the
    target depth. This should be one of three options: (1) `"unfilt"` will treat
    the sequencing depth of all (unfiltered) reads and sites as the 'full'
    depth; (2) `"mapqbaseq"` will filter out reads that don't pass the
    configured mapping or base quality, then calculate depth across all sites as
    the 'full' depth, (3) `"sitefilt"` will filter reads just as `"mapqbaseq"`
    does, but will only calculate the 'full' depth from sites passing the sites
    filter. As the main goal of subsampling is to make depth uniform for
    analyses, this latter option is preferred, as it will most accurately bring
    the depth of the samples to the target depth for analyses.
    (`"unfilt"`/`"mapqbaseq"`/`"sitefilt"`)
52. If `subsample_by` is set to `"unfilt"` or `"mapqbaseq"`,
    then it would be possible to recaculate extreme depth filters for the
    subsampled dataset. Enable this to do so, otherwise, the depth filters from
    the full depth bams will be used. If `subsample_by` is set to `"sitefilt"`
    this will have no effect, as the subsampling is already in reference to a
    set site list. (`true`/`false`)
53. When performing depth subsampling, you may want to leave some
    samples out that you kept in your 'full' dataset. These can be listed here
    and they will be removed from ALL depth subsampled analyses. A use case for
    this might be if you have a couple samples that are below your targeted
    subsample depth, and you don't want to include them. Note that if you
    configure multiple `subsample_dp`, these samples will be dropped from all of
    them. If you need to perform mutliple depth subsamplings with different
    subsets of samples, its best to run each depth individually. Alternatively,
    a config file can be made for each subsampled depth, however you may run
    into issues of file locking blocking both from running at the same time.
    (list of strings: `[]`)
54. Individually enable analyses to be performed with the
    subsampled data. These have the same function as described in the
    [`analyses`]section above. Enabling here
    will only run the analysis for the subsampled data, if you want to run it
    for the full data as well, you need to enable it in the analyses section as
    well. (`true`/`false`)
55. By default, this workflow will perform all analyses requested in the above
    section on all sites that pass the filters set in the above section. These
    outputs will contain `allsites-filts` in the filename and in the report.
    However, many times, it is useful to perform an analysis on different
    subsets of sites, for instance, to compare results for genic vs. intergenic
    regions, neutral sites, exons vs. introns, etc. Here, users can set an
    arbitrary number of additional filters using BED files. For each BED file
    supplied, the contents will be intersected with the sites passing the
    filters set in the above section, and all analyses will be performed
    additionally using those sites.

    For instance, given a BED file containing putatively neutral sites, one could
    set the following:

    ```yaml
    filter_beds:
      neutral: "resources/neutral_sites.bed"
    ```

    In this case, for each requested analysis, in addition to the
    `allsites-filts` output, a `neutral-filts` (named after the key assigned to
    the BED file in `config.yaml`) output will also be generated, containing the
    results for sites within the specified BED file that passed any set filters.
    More than one BED file can be set, up to an arbitrary number:

    ```yaml
    filter_beds:
      neutral: "resources/neutral_sites.bed"
      intergenic: "resources/intergenic_sites.bed"
      introns: "resources/introns.bed"
    ```

56. It may also sometimes be desireable to skip analyses on `allsites-filts`,
    say if you are trying to only generate diversity estimates or generate SFS
    for a set of neutral sites you supply.

    To skip running any analyses for `allsites-filts` and only perform them for
    the BED files you supply, you can set `only_filter_beds: true` in the config
    file. This may also be useful in the event you have a set of already
    filtered sites, and want to run the workflow on those, ignoring any of the
    built in filter options by setting them to `false`.

57. Phred-scaled mapping quality filter. Reads below this threshold will
    be filtered out. (integer)
58. Phred-scaled base quality filter. Reads below this threshold will be
    filtered out. (integer)
59. Determines whether overlapping read pairs will
    be clipped in BAM files supplied by users. This is useful as many variant
    callers will account for overlapping reads in their processing, but ANGSD
    will double count overlapping reads. If BAMs were prepped without this in
    mind, it can be good to apply before running through ANGSD. However, it
    essentially creates a BAM file of nearly equal size for every sample, so
    it may be nice to leave it off by applying it to the BAMs you feed into
    the pipeline. (`true`/`false`)
60. Additional options to pass to fastp trimming. (string)
    [[docs]](https://github.com/OpenGene/fastp)
61. Minimum overlap to collapse historical reads. Default
    in fastp is 30. This effectively overrides the `--length_required` option
    if it is larger than that. (INT) [[docs]](https://github.com/OpenGene/fastp)
62. Additional options to pass to bwa aln for mapping of historical
    sample reads. (string) [[docs]](https://bio-bwa.sourceforge.net/bwa.shtml)
63. Additional options to pass to Picard MarkDuplicates.
    `--REMOVE_DUPLICATES true` is recommended. (string)
    [[docs]](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)
64. Enable to run damage profiler on modern samples in
    addition to historical (`true`/`false`)
65. Parameters for pileup mappability analysis, see
    [GenMap's documentation](https://github.com/cpockrandt/genmap/) for more
    details.
66. Length of k-mers to calculate mappability for. (integer)
67. Allowed mismatches per k-mer. (integer)
68. A threshold mappability score. Each site gets an average
    mappability score taken by averaging the mappability of all K-mers that
    would overlap it. A score of 1 means all K-mers are uniquely mappable,
    allowing for `E` mismatches. This is done via a custom script, and may
    eventually be replaced by the SNPable method, which is more common.
    (float, 0-1)
69. Parameters for excluding sites based on extreme high
    and/or low global depth. The final sites list will contain only sites that
    pass the filters for all categories requested (i.e the whole dataset
    and/or the depth categories set in `samples.tsv`).
70. Whether you will generate extreme thresholds as a multiple of
    the median global depth (`"median"`) or as percentiles of the
    global depth distribution (`"percentile"`)
71. The bounds of the depth cutoff, defined as a numeric list. For
    the median method, the values will be multiplied by the median of the
    distribution to set the thresholds (i.e. `[0.5,1.5]` would generate
    a lower threshold at 0.5\*median and an upper at 1.5\*median). For the
    percentile method, these define the lower and upper percentiles to filter
    out (i.e [0.01,0.99] would remove the lower and upper 1% of the depth
    distributions). (`[FLOAT, FLOAT]`)
72. Whether to perform this filter on the dataset as a
    whole (may want to set to false if your dataset global depth distribution
    is multi-modal). (`true`/`false`)
73. Whether to perform this filter on the depth
    classes defined in the `samples.tsv` file. This will generate a global
    depth distribution for samples in the same category, and perform the
    filtering on these distributions independently. Then, the sites that pass
    for all the classes will be included. (`true`/`false`)
74. Seed to use when subsampling bams to lower depth.
    `"$RANDOM"` can be used to set a random seed, or any integer can be used
    to set a consistent seed. (string or int)
    [[docs]](https://www.htslib.org/doc/samtools-view.html)
75. General options in ANGSD, relevant doc pages are linked
76. Genotype likelihood model to use in calculation
    (`-GL` option in ANGSD, [[docs]](http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods))
77. When calculating individual depth, sites with depth higher
    than this will be binned to this value. Should be fine for most to leave
    at `100`. (integer, [[docs]](http://www.popgen.dk/angsd/index.php/Depth))
78. Individuals with sequencing depth below this value at a
    position will be treated as having no data at that position by ANGSD.
    ANGSD defaults to 1 for this. Note that this can be separately set for
    individual heterozygosity estimates with `mindepthind_heterozygosity`
    below. (integer, `-setMinDepthInd` option in ANGSD (INT)
    [[docs]](https://popgen.dk/angsd/index.php/Filters))
79. Used to fill the `-minInd` option for any dataset wide
    ANGSD outputs (like Beagles for PCA/Admix). Should be a floating point
    value between 0 and 1 describing what proportion of the dataset must have
    data at a site to include it in the output. (FLOAT)
    [[docs]](https://popgen.dk/angsd/index.php/Filters))
80. Used to fill the `-minInd` option for any population level
    ANGSD outputs (like SAFs or Beagles for ngsF-HMM). Should be a floating
    point value between 0 and 1 describing what proportion of the population
    must have data at a site to include it in the output. (FLOAT)
    [[docs]](https://popgen.dk/angsd/index.php/Filters))
81. Removes transitions using ANGSD, effectively removing them
    from downstream analyses. This is useful for removing DNA damage from
    analyses, and will automatically set the appropriate ANGSD flags (i.e.
    using
    [`-noTrans 1`](https://www.popgen.dk/angsd/index.php/SFS_Estimation#Brief_Overview)
    for SAF files and
    [`-rmTrans 1`](https://popgen.dk/angsd/index.php/Major_Minor) for Beagle
    files.) (`true`/`false`)
82. Additional options to pass to ANGSD during genotype likelihood
    calculation at all times. This is primarily useful for adding BAM input
    filters. Note that `--remove_bads` and `-only_proper_pairs` are enabled
    by default, so they only need to be included if you want to turn them
    off or explicitly show they are enabled. I've also found that for some
    datasets, `-C 50` and `-baq 1` can create a strong relationship between
    sample depth and detected diversity, effectively removing the benefits of
    ANGSD for low/variable depth data. I recommend that these aren't included
    unless you know you need them. Since the workflow uses BWA to map,
    `-uniqueOnly 1` doesn't do anything if your minimum mapping quality is
    \> 0, but you may wish to add it if you're bringing in your own files from
    another mapper. Mapping and base quality thresholds are also not needed,
    it will use the ones defined above automatically. If you prefer to correct
    for historical damage by trimming the ends of reads, this is where you'd
    want to put `-trim INT`. (string)
    (string, [[docs]](http://www.popgen.dk/angsd/index.php/Input#BAM.2FCRAM))
83. Same as `extra`, but only used when making SAF files (used
    for SFS, thetas, _F_~ST~, IBSrelate (not NgsRelate version),
    heterozygosity). Doesn't require options already in `extra` or defined via
    other params in the YAML (such as `notrans`, `minind`, `GL`, etc.)
    (string)
84. Same as `extra`, but only used when making Beagle and Maf
    files (used for PCA, Admix, ngsF-HMM, doIBS, NgsRelate, allele
    frequencies). Doesn't require options already in `extra` or defined via
    other params in the YAML (such as `rmtrans`, `minind`, `GL`, etc.)
    (string)
85. Method for inferring the major and minor alleles. Set to
    1 to infer from the genotype likelihoods, see
    [documentation](https://www.popgen.dk/angsd/index.php/Major_Minor)
    for other options. `1`, `2`, and `4` can be set without any additional
    configuration. `5` must also have an ancestral reference provided in the
    config, otherwise it will be the same as `4`. `3` is currently not
    possible, and is used for generating the dataset polarized MAFs. If you
    have a use for `3` that PopGLen is suited for, please open an issue and I
    will look into it. (integer)
86. The p-value to use for calling SNPs (float or string)
    [[docs]](http://www.popgen.dk/angsd/index.php/SNP_calling)
87. Method for inferring minor allele frequencies. Set to `1` to
    infer from genotype likelihoods using a known major and minor from the
    `domajorminor` setting above. Set to `2` to assume the minor is unknown.
    See [[docs]](http://www.popgen.dk/angsd/index.php/Allele_Frequencies) for
    more information. `3` is possible, and will estimate the frequencies both
    assuming a known and unknown minor. If you choose this option, you'll get
    both in the MAF outputs, but only the known will be passed to NgsRelate if
    you also use that. Other values are currently unsupported in PopGLen.
    (integer)
88. The minimum minor allele frequency required to call a SNP.
    This is set when generating the Beagle file, so will filter SNPs for
    PCAngsd, NGSadmix, ngsF-HMM, and NgsRelate. If you would like each tool
    to handle filtering for maf on its own you can set this to `-1`
    (disabled). (float, [[docs]](http://www.popgen.dk/angsd/index.php/Allele_Frequencies))
89. When estimating individual heterozygosity,
    sites with sequencing depth lower than this value will be dropped.
    (integer, `-setMinDepthInd` option in ANGSD) (integer)
90. Settings for ngsLD [[docs]](https://github.com/fgvieira/ngsLD)
91. For the LD estimates generated when setting
    `estimate_ld: true` above, set the maximum distance between sites in kb
    that LD will be estimated for (`--max_kb_dist` in ngsLD, integer)
92. For the LD estimates generated when setting
    `estimate_ld: true` above, randomly sample this proportion of pairwise
    linkage estimates rather than estimating all (`--rnd_sample` in ngsLD,
    float)
93. The same as `max_kb_dist_est-ld:`, but used when
    estimating LD decay when setting `ld_decay: true` above (integer)
94. The same as `rnd_sample_est-ld:`, but used when
    estimating LD decay when setting `ld_decay: true` above (float)
95. Additional plotting arguments to pass to
    `fit_LDdecay.R` when estimating LD decay (string)
96. When estimating LD decay, should the sample
    size corrected r^2^ model be used? (`true`/`false`, `true` is the
    equivalent of passing a sample size to `fit_LDdecay.R` in ngsLD using
    `--n_ind`)
97. The same as `max_kb_dist_est-ld:`, but
    used when linkage pruning SNPs as inputs for PCAngsd and NGSadmix. Pruning
    is performed on the whole dataset. Any positions above this distance will
    be assumed to be in linkage equilibrium during the pruning process.
    (integer)
98. The minimum r^2^ to assume two positions are
    in linkage disequilibrium when pruning for PCAngsd and NGSadmix analyses.
    (float)
99. Settings for NGSrelate
100. Any extra command line parameters to be passed to
     NgsRelate when performing an IBSrelate only (no allele frequencies) run.
     (string)
101. Any extra command line parameters to be passed to
     NgsRelate when performing a standard, allele frequency based, run.
     (string)
102. Settings for ngsF-HMM [[docs]](https://github.com/fgvieira/ngsF-HMM)
103. Set to `true` to run ngsF-HMM separately for each
     population in your dataset. Set to `false` to run for whole dataset at
     once. ngsF-HMM assumes Hardy-Weinberg Equilibrium (aside from inbreeding)
     in the input data, so select the option that most reflects this in your
     data. (`true`/`false`)
104. Whether or not to prune SNPs for LD before running the analysis.
     ngsF-HMM assumes independent sites, so it is preferred to set this to
     `true` to satisfy that expectation. (`true`/`false`)
105. The maximum distance between sites in kb
     that will be treated as in LD when pruning for the ngsF-HMM input. (INT)
106. The minimum r^2^ to assume two positions are in
     linkage disequilibrium when pruning for the ngsF-HMM input. Note, that
     this likely will be substantially higher for individual populations than
     for the whole dataset, as background LD should be higher when no
     substructure is present. (float)
107. Minimum ROH size in base pairs to include in inbreeding
     coefficient calculation. Set if short ROH might be considered low
     confidence for your data. (integer)
108. A list of integers that describe the size classes in base
     pairs you would like to partition the inbreeding coefficient by. This can
     help visualize how much of the coefficient comes from ROH of certain size
     classes (and thus, ages). List should be in ascending order and the first
     entry should be greater than `min_roh_length`. The first bin will group
     ROH between `min_roh_length` and the first entry, subsequent bins will
     group ROH with sizes between adjacent entries in the list, and the final
     bin will group all ROH larger than the final entry in the list. (list `[]`
     of integers)
109. Settings for realSFS
110. Whether or not to fold the produced SFS. Set to 1 if you have not
     provided an ancestral-state reference (0 or 1,
     [[docs]](http://www.popgen.dk/angsd/index.php/SFS_Estimation))
111. Determines number of bootstrap replicates to use when
     requesting bootstrapped SFS. Is used for both 1dsfs and 2dsfs (this is
     very easy to separate, open an issue if desired). Automatically used
     for heterozygosity analysis to calculate confidence intervals. (integer)
112. Settings for _F_~ST~ calculation in ANGSD
113. Determines which _F_~ST~ estimator is used by ANGSD. With 0
     being the default Reynolds 1983 and 1 being the Bhatia-Hudson 2013
     estimator. The latter is preferable for small or uneven sample sizes
     (0 or 1, [[docs]](http://www.popgen.dk/angsd/index.php/Fst))
114. Window size in bp for sliding window analysis (integer)
115. Window step size in bp for sliding window analysis (integer)
116. Settings for pi, theta, and Tajima's D estimation
117. Window size in bp for sliding window analysis (integer)
118. Window step size in bp for sliding window analysis (integer)
119. Minimum sites to include window in report plot. This does not
     remove them from the actual output, just the report plot and averages.
120. Settings for admixture analysis with NGSadmix. This analysis is
     performed for a set of K groupings, and each K has several replicates
     performed. Replicates will continue until a set of N highest likelihood
     replicates converge, or the number of replicates reaches an upper threshold
     set here. Defaults for `reps`, `minreps`, `thresh`, and `conv` can be left
     as default for most. Based on method described by
     [Pečnerová et al. 2021, _Curr. Biol._](https://doi.org/10.1016/j.cub.2021.01.064)
121. A list of values of K to fit the data to (list `[]` of
     integers)
122. The maximum number of replicates to perform per K. Default is 100.
     (integer)
123. The minimum number of replicates to perform, even if
     replicates have converged. Default is 20. (integer)
124. The convergence threshold - the top replicates must all be
     within this value of log-likelihood units to consider the run converged.
     Default is 2. (integer)
125. The number of top replicates to include in convergence
     assessment. Default is 3. (integer)
126. Additional arguments to pass to NGSadmix (for instance,
     increasing `-maxiter`). (string,
     [[docs]](http://www.popgen.dk/software/index.php/NgsAdmix))
127. Settings for identity by state calculation with ANGSD
128. Whether to use a random (1) or consensus (2) base in IBS
     distance calculation
     ([[docs]](http://www.popgen.dk/angsd/index.php/PCA_MDS))
