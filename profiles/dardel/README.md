# Snakemake profile for running PopGLen on PDC's Dardel system

PopGLen was developed on UPPMAX's Rackham and PDC's Dardel systems. As Rackham
is now end of life, it continues to be developed on Dardel. Here is the profile
I have used to run the workflow on Dardel, with some of my preferred settings.

This profile is intended to be passed to Snakemake when the pipeline is run
using the `--profile` option. Since it is submitting jobs to the SLURM queue, it
is intended to be run in a screen with few local cores, so as to not negatively
impact the login node.

## What is specific to Dardel in here?

This will automatically reserve 1700MB of memory per thread, the standard on
Dardel. Rules that are likely to run out of memory are set to rerun with more
threads if they fail, so this should scale well. You may still need to increase
resources for certain rules depending on your data.

The SLURM partition is set to 'shared' the standard for jobs that use partial
nodes like ours will. We also set nodes to 1 so that we ensure the threads are
all reserved on one node (the memory reservation should do this anyway, so this
is just to be safe).

We by default pass `--tmp-sandbox -B /cfs/klemming` to any Singularity commands,
which is needed for the containers to work properly on Dardel.

> [!WARNING]
>
> **Use Singularity with this profile, not Apptainer**
>
> Apptainer no longer has `--tmp-sandbox` as an option for its `exec` command.
> As such, these arguments won't work if you are using Apptainer rather than
> Singularity. For that reason, if you're using this profile on Dardel, be
> sure to either install singularity in your conda environment or load the
> Singularity module: `ml load PDC singularity`

**What you need to change:**

- `conda-prefix:` Either remove this, or set it to a folder you have created
  where your conda environments for Snakemake are going to be stored. This is
  useful if you will run PopGLen multiple times in different working dirs, as
  the environments will be shared across them, saving space and creation time on
  initial startup. I find this is a good way to use the storage on your compute
  project as the environments don't take up much space, but use a lot of files,
  which the compute storages are tailored for. If you remove this line, the
  environments will be stored in the working directory in the `.snakemake`
  folder.
- `singularity-prefix:` The same as for `conda-prefix`, but for singularity
  containers. I recommend setting it in a folder next to the conda one.
- `slurm_account=` Change this to the compute account you want your compute
  hours to be billed to.
- `tmpdir=` Change this to your user's allocated temporary directory. This
  should be at `/cfs/klemming/scratch/first-letter-of-username/username`.
