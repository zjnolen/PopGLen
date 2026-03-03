import pandas as pd
import sys

# Python script to concatenate multiple files and optionally add a header.
# `inputs` is a list of files to concatenate (will be treated as tab-separated
# files without header), `output` is the concatenated output file (will be made
# as tab separated with header), and `header` is a list of column names to add
# as a header (or set to false to not add a header to the output, set to false
# automatically if the rule calling this script has no header param).

sys.stderr = open(snakemake.log[0], "w") if snakemake.log else sys.stderr

inputs = snakemake.input
output = snakemake.output[0]

# Check if snakemake rule has a header to use in the params, if not, no header
# will be output, even if it exists in the input files.
if hasattr(snakemake.params, "header"):
    header = snakemake.params.header
else:
    header = False

# Check if snakemake rule specifies input files already have header. If so, skip
# the first line of each input.
if hasattr(snakemake.params, "input_has_header") and snakemake.params.input_has_header:
    skiprows = 1
else:
    skiprows = None

df = pd.concat(
    [
        pd.read_csv(f, sep="\t", dtype=str, header=None, skiprows=skiprows)
        for f in inputs
    ],
    ignore_index=True,
)

df.to_csv(output, sep="\t", index=False, header=header)
