# Profiles

These are some example profiles you can use with the pipeline. Generally,
two profiles get applied when using Snakemake, the 'global profile' and the
'workflow-profile'. Read more about the difference between global and workflow
specific profiles in Snakemake
[here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).

The 'global profile' (`--profile`) is mostly for general options to make
Snakemake run on your system, things like your machine's core count or your
slurm deployment. Take a look at the dardel profile here for an example of
one I use on PDC's Dardel.

The 'workflow-profile' (`--workflow-profile`) contains workflow-specific options,
usually modifying rule resources. A profile called 'default' is here and will be
applied when running the workflow by default. This means **if you want to change
resources, you must set them in the command line, edit `profiles/default/config.yaml`,
or pass them in a profile to `--workflow-profile`**. If you add them in the profile
you set with `--profile` or change them in the `.smk` files, they will end up overwritten
by the default workflow profile.
