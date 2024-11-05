# Running PopGLen on a cluster using a job queue

PopGLen is primarily made for running on a high performance computing system
that has a job scheduler. This system allows for each rule to be submitted as a
single job on a system with a large number of nodes, enabling large levels of
parallelization. Snakemake is well integrated to many job scheduling systems,
and some examples can be found on the
[Snakemake plugin catalog](https://snakemake.github.io/snakemake-plugin-catalog/index.html).

## Command line options for cluster execution

Here, we will walk through what an execution might look like on an HPC with the
SLURM job scheduler. This requires a conda environment with Snakemake installed,
along with the snakemake-executor-plugin-slurm installed:

```bash
conda create \
    -n popglen \
    -c conda-forge -c bioconda \
    snakemake=8.25.0 \
    snakemake-executor-plugin-slurm

conda activate popglen

## OR, if you already set up a popglen environment and just need to add the
## executor plugin:

conda activate popglen

conda install -c conda-forge -c bioconda snakemake-executor-plugin-slurm
```

Once you've set up a working directory and workflow for PopGLen, you can run it
using the following command:

```bash
snakemake --use-conda --use-singularity --executor slurm \
    --set-resources slurm_account=<project> slurm_partition=<partition name>
```

This will ensure that rules will be submitted as jobs to the job queue, by
default using the SLURM account and partition defined in `--set-resources`. The
values entered for these will depend on your system, but would be the equivalent
of the account and partition options you set in the header of SLURM scripts:

```bash
#SBATCH -A <project>
#SBATCH -p <partition name>
```

## Combining command line options into a profile

You may wish to set a few more options, such as the maximum number of threads
available to each job, the maximum number of jobs to have running and in the
queue at once based on your queue's limits, and how many local cores to use for
rules that run outside the job system (in PopGLen these are quick commands like
making lists of BAM files or making symbolic links). These quickly begin to add
up, so take a look at how we define these all in an example profile for the HPC
system Dardel at PDC:

```yaml title="profiles/dardel/config.yaml"
restart-times: 3
local-cores: 1
use-conda: true
use-singularity: true
jobs: 999
keep-going: true
max-threads: 128
executor: slurm
singularity-args: '--tmp-sandbox -B /cfs/klemming'
default-resources:
  - "mem_mb=(threads*1700)"
  - "runtime=60"
  - "slurm_account=<account>"
  - "slurm_partition=shared"
  - "tmpdir='/cfs/klemming/scratch/u/user'"
```

Assuming we place this file in the working directory under `profiles/dardel`,
we can run Snakemake with simply:

```bash
snakemake --profile ./profiles/dardel
```

and it will automatically use the options defined in the profile. Here, we allow
up to 3 retries for failed jobs, 1 local core for running local jobs, enable
both conda and singularity usage, set a max of 999 jobs to be running at once,
tell rules to keep going if one fails, set a maximum of 128 threads for a single
job (the size of a node on Dardel), define SLURM as the executor, set some
required arguments to be passed to Singularity on Dardel, and set several
default resources, including giving each rule 1.7GB memory per thread, which is
the amount available on Dardel's shared partition. We also set the temporary
directory to use, as it is not set on Dardel by default.

Note, that if you're not on Dardel, this will need some changes, so make a new
profile that matches your system. You'll likely need to change `max-threads`,
`slurm_account`, `slurm_partition`, and `tmpdir`. `mem_mb` can have 1700 changed
to match the amount of memory your system has per core. `singularity-args` is
specific to Dardel, and can be omitted, unless your system requires something
similar. At the minimum, `-B /cfs/klemming` will need to be changed or removed.

## Executing the workflow while away using Screen

If you run Snakemake in the login shell of your system, it will cancel when you
logout, which is not ideal for the expected long runtimes for processing WGS
data. If your system lets you submit new jobs from within other jobs, you can
submit your Snakemake command as a long running SLURM job. Be sure to activate
the conda environment for PopGLen either inside the job before running Snakemake
or before submitting the job, as the job inherits your active environment.

If you can't submit jobs from inside other jobs, or you would like the option to
interact with it while its running, you can run it inside a screen, if screen is
available on your system. First start a new screen:

```bash
screen -S project-name
```

This will open up a virtual terminal that will stay active as long as your
system is online. Inside this terminal, you can activate snakemake, do dry runs,
and start the workflow:

```bash
conda activate popglen

# do a dry run
snakemake --profile ./profiles/my-cluster -n

# do a real run
snakemake --profile ./profiles/my-cluster
```

Snakemake will then submit jobs from inside the screen. You can disconnect from
the screen with CTRL-A + D, and safely log out without Snakemake being
interrupted. Then, you can go back in by resuming your screen:

```bash
screen -r project-name
```

and see how far it has gotten (or if it fails and you need to change something).

When you're done with a screen, you can kill it with CTRL-A + K.
