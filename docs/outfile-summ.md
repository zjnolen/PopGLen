# Output files

Here is a brief overview of the file paths for output files. Files marked with
*temp* are flagged as temporary and will be deleted when they are no longer
needed by the workflow. Paths have wildcards in their names, indicated by `{}`.
Some of these wildcards include:

- `{dataset}` - The dataset name defined in the config file
- `{ref}` - The reference genome name defined in the config file
- `{sample}` - Sample names as written in the sample list
- `{population}` - Population names as written in the sample list, but also
  including depth groupings in some rare cases
- `{unit}` - The sequencing unit (usually flow cell barcode + sequencing lane)
  of the sequencing reads
- `{lib}` - The sequencing library the sequencing reads come from
- `{dp}` - The subsampled sequencing depth, if enabled. The full depth results
  will have this wildcard empty in the path.
- `{sites}` - The filter set used to produce the result. If no user provided
  filters are given, this will always be `allsites`. If user provided filters
  are in the config, it will use the name given to the BED file.

## Mapping

### Trimmed/Collapsed reads

#### Trimmed reads (uncollapsed, *temp*)

`results/preprocessing/fastp/{sample}_{unit}_{lib}.{read}.paired.fastq.gz`

#### Collapsed reads (*temp*)

`results/preprocessing/fastp/{sample}_{unit}_{lib}.merged.fastq.gz`

### Processed alignments

#### BAM files

`results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.bam` (also symlinked to
`results/datasets/{dataset}/bams/{sample}.{ref}.bam`)

#### BAM files (MapDamage rescaled)

`results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.rescaled.bam` (also
symlinked to `results/datasets/{dataset}/bams/{sample}.{ref}.bam` when enabled)

#### Subsampled BAM files

`results/datasets/{dataset}/bams/{sample}.{ref}{dp}.bam`

## Quality control

### fastp reports

`results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_merged.html` (collapsed
reads)

`results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_paired.html`
(trimmed, i.e. uncollapsed, reads)

### Picard MarkDuplicates metrics

`results/mapping/qc/mark_duplicates/{sample}.{ref}.picard.metrics`

### Dedup duplicate removal metrics

`results/mapping/qc/dedup/{sample}_{lib}.{ref}.dedup.json`

### Bamutil overlap clipping metrics

`results/mapping/qc/bamutil_clipoverlap/{sample}.{ref}.rmdup.realn.clipoverlap.stat`

### Samtools flagstat for raw alignments (pre-duplicate removal)

`results/mapping/mapped/{sample}_{library}.{ref}.merged.flagstat` (collapsed
reads)

`results/mapping/mapped/{sample}.{ref}.paired.flagstat` (trimmed, i.e.
uncollapsed, reads)

### Qualimap reports

`results/mapping/qc/qualimap/{sample}.{ref}/qualimapReport.html` (BAMs
processed by PopGLen),

`results/datasets/{dataset}/qc/user-provided-bams/qualimap/{sample}.{ref}/qualimapReport.html`
(BAMs brought by user),

`results/datasets/{dataset}/qc/qualimap/qualimap_all.{ref}_mqc.html` (MultiQC)

### Mean and st. dev. depth, mapping rates

`results/datasets/{dataset}/qc/{dataset}.{ref}_all{dp}.sampleqc.tsv`

### Identity by state distance to reference

`results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_{filts}.refibs.tsv`
(here `{filts}` refers to either unfiltered, or with one of the sites filters)

`results/datasets/{dataset}/plots/ibs_refbias/{dataset}.{ref}_all{dp}_{filts}.{grouping}.pdf` (plots)

### DNA damage

`results/mapping/qc/damageprofiler/{sample}.{ref}` (DamageProfiler)

`results/mapping/qc/mapdamage/{sample}.{ref}` (MapDamage 2)

`results/datasets/{dataset}/qc/dna-damage-mqc/dna-damage_all.{ref}_mqc.html`
(MultiQC)

## Filtering

### Reference links and indices

`results/ref/{ref}/{ref}.fa` (symlink to provided referenceindices are built off
this symlink)

`results/ref/{ref}/{ref}.anc.fa` (symlink to provided ancestral state reference,
if provided in config)

### Reference 'chunk' region files

`results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf` (where
`{chunk}` is the chunk number)

#### Mappability

##### Genmap mappability scores

`results/ref/{ref}/genmap/map/{ref}_k{k}_e{e}.bedgraph`

##### Pileup mappability scores

`results/ref/{ref}/genmap/pileup/{ref}_pileup_mappability_k{k}_e{e}.bed`

##### Pileup mappability filter BED

`results/datasets/{dataset}/filters/pileupmap/{ref}_k{k}_e{e}_{thresh}.bed`

#### Repeats

##### Repeat library generated in pipeline

`results/ref/{ref}/repeatmodeler/{ref}-families.fa`

##### Identified repeats (either from provided or generated library)

`results/ref/{ref}/repeatmasker/{ref}.fa.out.gff`

##### Repeat filter BED

`results/ref/{ref}/repeatmasker/{ref}.fa.out.bed`

#### Combined filter BED

`results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_allsites-filts.bed`

`results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.bed`
(when user provides additional filter BEDs. `{sites}` will be the name in the
config of the filter BED)

## Genotype likelihood analyses

### GL files

#### Beagles

`results/datasets/{dataset}/beagles/{dataset}.{ref}_{population}{dp}_{sites}-filts.beagle.gz`

#### SAFs

`results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf.idx`

`results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf.pos.gz`

`results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf.gz`

(for the above files, `{population}` can also be from samples for 1 sample SAF)

### Linkage disequilibrium

#### LD estimates

`results/datasets/{dataset}/analyses/ngsLD/{dataset}.{ref}_{population}{dp}_{sites}-filts.ld_maxkbdist-{maxkb}_rndsample-{rndsmp}.gz`
(where `{maxkb}` is the max distance between SNPs in kilobases and `{rndsnp}` is
the proportion of SNPs to sample for the LD calculation)

#### Relatedness tables

##### IBSrelate (SFS-based)

`results/datasets/{dataset}/analyses/kinship/ibsrelate_sfs/{dataset}.{ref}_all{dp}_{sites}-filts.kinship`

##### IBSrelate (IBS-based)

`results/datasets/{dataset}/analyses/kinship/ibsrelate_ibs/{dataset}.{ref}_all{dp}_{sites}-filts.kinship`

##### IBSrelate (NgsRelate implementation)

`results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts_ibsrelate-nofreq.tsv`

##### NgsRelate (allele frequency-based)

`results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts_ngsrelate-freq.tsv`

#### Principal component analysis

##### Covariance matrix

(for these, `{population}` is either `all` or `all_excl_pca-admix` depending on
if relatives were removed)

`results/datasets/{dataset}/analyses/pcangsd/{dataset}.{ref}_{population}{dp}_{sites}-filts.cov`

#### Admixture

##### NGSadmix outputs

(for these, `{population}` is either `all` or `all_excl_pca-admix` depending on
if relatives were removed)

`results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.qopt`

`results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.fopt.gz`

`results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}_optimization_wrapper.log`
(this logs the replicate number, seed, and log-likelihood of the replicate
NGSadmix runs for each value of K)

##### EvalAdmix correlation of residuals

(for these, `{population}` is either `all` or `all_excl_pca-admix` depending on
if relatives were removed)

`results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_{population}{dp}_{sites}-filts_K{kvalue}.corres`

#### Site frequency spectrum (SFS)

##### One population/one sample

`results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.sfs`

`results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.boot.sfs`
(bootstrapped)

##### Two population/two sample

`results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.sfs`

`results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.boot.sfs`
(bootstrapped)

#### Theta and neutrality statistics (nucleotide diversity, Watterson's theta, Tajima's D)

##### Windowed theta estimates

`results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_{population}{dp}_{sites}-filts.thetaWindows.{win}_{step}.pestPG`

##### Mean and confidence intervals

`results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_all{dp}_{sites}-filts.window_{win}_{step}.{stat}.mean.tsv`
(where `{stat}` is `pi`, `watterson`, or `tajima`)

#### Genetic differentiation (*F*~ST~)

##### Windowed *F*~ST~ estimates

`results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.fst.window_{win}_{step}.tsv`

##### Global estimates

`results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{unit}pairs{dp}_{sites}-filts.fst.global.tsv`
(where `{unit}` is either `pop` or `ind`, depending on whether *F*~ST~ was
requested between populations or individuals)

#### Individual heterozygosity table

`results/datasets/{dataset}/analyses/heterozygosity/{dataset}.{ref}_all{dp}_{sites}-filts_heterozygosity.tsv`

#### Inbreeding

##### ngsF-HMM model based inbreeding coefficients

`results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts.indF`

##### IBD segments (runs of homozygosity)

`results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_{population}{dp}_{sites}-filts.roh`
(raw)

`results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_all{dp}_{sites}-filts.all_roh.bed`
(after filtering out runs less than minimum length in configuration file)

##### *F~RoH~* table

`results/datasets/{dataset}/analyses/ngsF-HMM/{dataset}.{ref}_all{dp}_{sites}-filts.ind_froh.tsv`

#### Identity by state matrix

`results/datasets/{dataset}/analyses/IBS/{dataset}.{ref}_all{dp}_{sites}-filts.ibsMat`

#### Population specific allele frequencies

`results/datasets/{dataset}/mafs/{dataset}.{ref}_{population}{dp}_{sites}-filts.pop-maj.mafs.gz`
(sites segregating in population, major/minor polarized to population freqs
unless using `-doMajorMinor 4` or `5`)

`results/datasets/{dataset}/mafs/{dataset}.{ref}_{population}{dp}_{sites}-filts.dataset-maj.mafs.gz`
(sites segregating in dataset, including invariant in population, major/minor
polarized to dataset wide frequencies unless using `-doMajorMinor 4` or `5`)
