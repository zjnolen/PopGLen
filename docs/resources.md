# Altering resource usage

This workflow sets a default number of threads and runtime for each rule. These
defaults are set both in the rules themselves and given in the default workflow
profile found in `profiles/default/config.yaml`. Memory is not set per rule as
development took place on clusters that did not allow direct requests of memory
when submitting SLURM jobs, instead always allocating 6.4GB RAM per core for each
job. The easiest way to get up and running with this workflow is to set
`--default-resources mem_mb="XXXX * threads"` when running the Snakemake command,
replacing `XXXX` with the RAM available per core on your HPC system in MB (in
the case of ours, this was 6400). See the
[profiles](https://github.com/zjnolen/PopGLen/tree/v0.4.3/profiles) in the GitHub
repository as an example.

The default resources may not always work, your data may need more
memory, or longer runtimes, or maybe you even need shorter if your HPC has
shorter runtime limits than some of the defaults we set (up to 3 days).
Snakemake makes it easy to alter resources in the command line using the
`--set-resources` and `--set-threads` options, which will override anything set
in the workflow already.

Instead of using command line arguments to change the resources, you can also
edit the default workflow profile or replace it using `--workflow-profile`
(in addition to any global profiles you set with `--profile`). You can find the
default profile to edit or use as a template in
[profiles/default](https://github.com/zjnolen/PopGLen/tree/v0.4.3/profiles/default).
It is important to remember that for Snakemake, command line arguments take precedence
over the workflow profile (`--workflow-profile`), which takes precedence over the
global profile (`--profile`). So, if you plan to edit resources, do it in the command
line or a `--workflow-profile`, as the default workflow profile will overwrite any 
resources set in `--profile`.

Here is an example of what it would look like to change the number of threads
for the rule `bwa_index` from the default of 1 to a new value of 5:

```yaml linenums="1" hl_lines="5" title="profiles/default/config.yaml"
set-threads:
  # Reference Prep
  link_ref: 1
  link_anc_ref: 1
  bwa_index: 5
  samtools_faidx: 1
  ref_chunking: 1
  picard_dict: 1
```

And if you wanted to update the runtime to give it only a maximum of 1 day and
add a limit of 2GB for memory:

```yaml linenums="150" hl_lines="7-9" title="profiles/default/config.yaml"
set-resources:
  # Reference Prep
  link_ref:
    runtime: "5m"
  link_anc_ref:
    runtime: "5m"
  bwa_index:
    runtime: "1d"
    mem_mb: 2048
  samtools_faidx:
    runtime: "1h"
  ref_chunking:
    runtime: "5m"
  picard_dict:
    runtime: "10m"
```

You can also use the term 'attempt' in these definitions, which allow you to
scale the resources with the number of attempts a rule has made, automatically
increasing threads, runtime, or memory with each attempt.
