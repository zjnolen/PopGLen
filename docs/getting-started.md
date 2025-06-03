# Getting Started

## Tutorial

[A tutorial is available](tutorial.md) with a small(ish) dataset where
biologically meaningful results can be produced. This can help get an
understanding of a good workflow to use different modules. You can also follow
along with your own data and just skip analyses you don't want. I recommend at
least looking at this to get an idea of the scope of the pipeline as well as to
see a few tips. If you prefer to just jump in instead, below describes how to
get a new project up and running.

## Requirements

This pipeline can be run on Linux systems with Conda and Apptainer/Singularity
installed. All other dependencies will be handled with the workflow, and thus,
sufficient storage space is needed for these installations (~10GB, though this
is only when using *all* features of the pipeline, if you are not using any of
the mapping module, it will be more like 3-5GB). It can be run on a local
workstation with sufficient resources and storage space (dataset dependent), but
is aimed at execution on high performance computing systems with job queuing
systems.

Data-wise, you'll need a reference genome (uncompressed) and some sequencing
data for your samples. The latter can be either raw FASTQ files, BAM alignments
to the reference, or SRA accession numbers for already published FASTQ files.

## Deploying the workflow

The pipeline can be deployed in two ways: (1) using
[Snakedeploy](https://github.com/snakemake/snakedeploy) which will deploy the
pipeline as a module (recommended); (2) clone the repository at the
version/branch you prefer (recommended if you will need to change workflow code
beyond what is possible with module definitions). If you are curious how
modularization works in Snakemake, take a look at the
[docs](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).

Both methods require a Snakemake environment to run the pipeline in.

### Preparing the environment

First, create an environment for Snakemake, including Snakedeploy if you intend
to deploy that way:

```bash
conda create -c conda-forge -c bioconda --name popglen snakemake snakedeploy
```

If you already have a Snakemake environment, you can use that instead. Snakemake
versions >=8 are likely to work, but most testing is on 8.20. If you intend to
use a job queuing system with the pipeline, be sure to install the appropriate
[executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/). Most
of the testing has been with `snakemake-executor-plugin-slurm`, which can be
installed in the environment you are running Snakemake from.

Activate the Snakemake environment:

```bash
conda activate snakemake
```

### Option 1. Deploying with Snakedeploy

!!! warning "Warning for clusters where worker nodes have no network access."
    If you will use a job queue on a cluster where worker nodes have no network
    access, Snakedeploy will not work properly, as rules with external scripts
    retrieve their code inside the submitted job, requiring network access when
    deployed as a module. If this is your configuration, you should clone the
    repository instead (see Option 2 below). Note this will also affect rules in
    Snakemake wrappers, so see the [cluster configuration docs](cluster.md) for
    more information on usage on clusters with this configuration.

Make your working directory:

```bash
mkdir -p /path/to/work-dir
cd /path/to/work-dir
```

And deploy the workflow, using the tag for the version you want to deploy:

```bash
snakedeploy deploy-workflow https://github.com/zjnolen/PopGLen . --tag v0.4.1
```

This will generate a simple Snakefile in a `workflow` folder that loads the
pipeline as a module. It will also download the template `config.yaml`,
`samples.tsv`, and `units.tsv` in the `config` folder.

### Option 2. Cloning from GitHub

Go to the folder you would like you working directory to be created in and
clone the GitHub repo:

```bash
git clone https://github.com/zjnolen/PopGLen.git
```

If you would like, you can change the name of the directory:

```bash
mv PopGLen work-dir-name
```

Move into the working directory (`PopGLen` or `work-dir-name` if you changed it)
and checkout the version you would like to use:

```bash
git checkout v0.4.1
```

This can also be used to checkout specific branches or commits.

## Configuring the workflow

Now you are ready to configure the workflow, see the documentation for that
[in the configuration section](config.md).
