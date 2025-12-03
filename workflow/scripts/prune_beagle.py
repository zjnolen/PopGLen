import sys
import pandas as pd
import gzip
from csv import DictReader

sys.stderr = open(snakemake.log[0], "w")


def prune_beagle(beagle, pos):
    """
    Subsets a beagle file (gzipped) to only include sites in the pos file.

    beagle - path to beagle file to be pruned (gzipped)
    pos - path to file with positions to keep (one per line, format: chr:pos)
    """
    with gzip.open(beagle, "rt") as f:
        beagle_df = pd.read_csv(f, sep="\t", dtype=str, header=None)
    beagle_df.columns = beagle_df.iloc[0]
    beagle_df = beagle_df.drop(0)
    with open(pos, "r") as f:
        pos_list = [i.replace(":", "_") for i in f.read().splitlines()]
    pruned_beagle_df = beagle_df[beagle_df.iloc[:, 0].isin(pos_list)]
    if len(pruned_beagle_df.index) != len(pos_list):
        raise ValueError(
            "Number of sites in pruned beagle file does not match the number of"
            "sites in the pos file."
        )
    return pruned_beagle_df


def df2gzip(df, output_path):
    """
    Writes a dataframe to a gzipped tsv file.
    """
    with gzip.open(output_path, "wt") as f:
        df.to_csv(f, sep="\t", index=False)


pruned_beagle_df = prune_beagle(snakemake.input.beagle, snakemake.input.sites)
df2gzip(pruned_beagle_df, snakemake.output.beagle)
