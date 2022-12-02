# General Settings

Here you will find example configuration files for this pipeline. For your 
run, you'll need to edit `samples.tsv`, `units.tsv`, and `config.yaml`.

## Samples list

All your samples should be listed in `samples.tsv`. Samples preceded with a 
`#` will not be included, this can be useful if you want to exclude a sample 
after quality checking.

Each sample must have four columns filled in the sample sheet. The columns are 
tab separated:

- `sample` - The name of the sample
- `population` - The populations you will group your samples into. These are 
  the groups that population level stats are calculated on.
- `time` - This should be either `modern` or `historical`, the only thing this 
  will affect is whether or not your bam files will be corrected for 
  post-mortem damage.
- `depth` - This is only used for depth filtering. Extreme levels (both high 
  and low) will be calculated for each group you put here as well as the 
  dataset as a whole, and all will be filtered out for all analyses.

The values in the sample list will end up in filenames, so it is best if they 
are characters that do not need escaping and do not sometimes run into issues 
with text encoding. I recommend only using alphanumeric characters. 
Non-English alphabet letters should work, especially on systems configured for 
them, but have not been fully tested.

## Units list

All your raw data will be pointed to in `units.tsv`.

Each sample must have four columns filled in the sample sheet. The columns are 
tab separated:

- `sample` - The sample name, same as in `samples.tsv`.
- `unit` - This is used to fill out the ID read group in the bam file. 
  Usually, it will be a sequencer barcode and land number in the format 
  `barcode.lane#`. This doesn't matter much if you have only sequenced your 
  samples one time each, and this pipeline doesn't support merging multiple 
  sequencing batches per samples, so it is safe to just assign anything to 
  this for now, its more for future improvements.
- `platform` - This is used to fill out the PL read group. Put what you'd want 
  there. Usually `ILLUMINA` for Illumina platforms.
- `fq1` and `fq2` - The absolute or relative paths from the working directory 
  to the raw fastq files for the sample. Currently the pipeline only supports 
  paired-end sequencing, single end may be added down the line.

## Configuration file

`config.yaml` contains the configuration for the workflow, this is where 
you will put what analyses, filters, and options you want. See the comments in 
the configuration file for information on each option.

### Recommended workflow configuration

There is a recommended order in which to do run the analyses so that you may 
check in at key points in the pipeline. This prevents running many downstream 
analyses before ensuring that earlier settings are as you want them.

For all runs, you'll want to fill out the dataset and reference configuration.

#### Run 1 - Sample processing, checking sample and filter quality

The first run is recommended to get all your samples mapped, quality estimates 
made, and to determine if the filters you've set are correct. For this, it is 
recommended your analyses section looks something like this:

```yaml
analyses:

# filtering
  # filter out regions with low mappability (true/false)
  genmap: true
  # filter out repetitive regions. fill in for only one of the three options following
  repeatmasker:
    # path to a local repeat library .fa file
    local_lib:
    # string with taxa name to get repeat libraries from dfam
    dfam_lib:
    # run repeatmodeler to build repeat database (true/false)
    build_lib: true
  # filter out sites at the extremes of the depth distribution (upper and lower 2.5%)
  extreme_depth: true
  # filter sites with data in less than this proportion of samples [0-1]
  dataset_missing_data: 0.9
  # filter sites with data in less than this proportion of samples per population [0-1]
  population_missing_data: 0.8

# quality control (all true/false)
  # generate qualimap report for all bam files
  qualimap: true
  # generate damage report for historical samples
  damageprofiler: true
  # generate endogenous content proportions for all samples
  endogenous_content: true
  # generate a table of relatedness between all individuals
  relatedness: true

# population genomic analyses (all true/false)
  # PCA
  pca_pcangsd: false
  # Admixture
  admix_ngsadmix: false
  # pi, wattersons, and tajima's d
  thetas_angsd: false
  # individual heterozygosity
  heterozygosity_angsd: false
  # pairwise population fst (only set `populations` to true/false)
  fst_angsd:
    populations: false
  # inbreeding using runs of homozygosity
  inbreeding_ngsf-hmm: false 
  # identity by state matrix
  ibs_matrix: false
```

Basically, you will have all your filters set up, and your quality checks 
turned on, but no analyses. Run the pipeline, generate a report, and see if 
things are matching what you'd expect - i.e. no samples are of excessively low 
quality, misidentified, highly related to others, and that your filters do not 
seem to be too loose or stringent. If you need to make adjustments, make them, 
rerun the pipeline, and settle on the values for the first section.

Note that the first run takes the longest, especially if you are using 
RepeatModeler. This is because mapping and RepeatModeling are by far the 
longest running parts of this pipeline.

#### Run 2 - Run your analyses

The next run can be for any combination of analyses. Switch whichever analyses 
you want to run to `true` and re-run, only what is needed to make these now 
will be run. Then you can generate another report with them included.