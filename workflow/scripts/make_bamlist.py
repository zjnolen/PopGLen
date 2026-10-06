import os
import sys

# Prepares a bamlist for use as input to ANGSD. This is a list of bam files, one
# per line given as the absolute path. There is no newline at the end of the
# file. These concerns may no longer be present in the current version of ANGSD,
# but are maintained in this script.

sys.stderr = open(snakemake.log[0], "w") if snakemake.log else sys.stderr

inputs = snakemake.input.bams
bamlist = snakemake.output[0]

print(f"BAM order: {snakemake.params.sampord}", file=sys.stderr)

abs_paths = "\n".join([os.path.realpath(b) for b in inputs])

with open(bamlist, "w") as f:
    f.write(abs_paths)
