import pandas as pd

# Python script to concatenate multiple files and optionally add a header.
# `inputs` is a list of files to concatenate (will be treated as tab-separated
# files without header), `output` is the concatenated output file (will be made
# as tab separated with header), and `header` is a list of column names to add
# as a header (or set to false to not add a header to the output, set to false
# automatically if the rule calling this script has no header param).

inputs = snakemake.input
output = snakemake.output[0]
if hasattr(snakemake.params, "header"):
    header = snakemake.params.header
else:
    header = False

df = pd.concat(
    [pd.read_csv(f, sep="\t", dtype=str, header=None) for f in inputs],
    ignore_index=True,
)

df.to_csv(output, sep="\t", index=False, header=header)
