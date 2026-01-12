# Altering resource usage

This workflow sets a default number of threads and runtime for each rule. Memory
is not set per rule as development took place on clusters that did not allow
direct requests of memory when submitting SLURM jobs, instead always allocating
6.4GB RAM per core for each job. The easiest way to get up and running with this
workflow is to set `--default-resources mem_mb="XXXX * threads"` when running
the Snakemake command, replacing `XXXX` with the RAM available per core on your
HPC system in MB (in the case of ours, this was 6400). See the
[profiles](https://github.com/zjnolen/PopGLen/tree/v0.4.3/profiles) in the
GitHub repository as an example.

However, the default resources may not always work, your data may need more
memory, or longer runtimes, or maybe you even need shorter if your HPC has
shorter runtime limits than some of the defaults we set (up to 7 days).
Snakemake makes it easy to alter resources in the command line using the
`--set-resources` and `--set-threads` options, which will override anything set
in the workflow already.

This is even better set up in a profile. PopGLen includes a default profile in
the
[profiles/default](https://github.com/zjnolen/PopGLen/tree/v0.4.3/profiles/default)
folder, which has all the rules already populated with the default resources
already set. This profile is downloaded whether you cloned the repository or
used Snakedeploy, so all you need to do is change the values to match your needs
(which may require adding an additional resource option, like `mem_mb`).

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

```yaml linenums="150" hl_lines="5-7" title="profiles/default/config.yaml"
set-resources:
#   # Reference Prep
  link_ref:
    runtime:
  link_anc_ref:
    runtime:
  bwa_index:
    runtime: "1d"
    mem_mb: 2048
  samtools_faidx:
    runtime:
  ref_chunking:
    runtime:
  picard_dict:
    runtime:
```

You can also use the term 'attempt' in these definitions, which allow you to
scale the resources with the number of attempts a rule has made, automatically
increasing threads, runtime, or memory with each attempt. This is already done
for several rules in the config, largely to automatically request more memory
due to OOM errors.
