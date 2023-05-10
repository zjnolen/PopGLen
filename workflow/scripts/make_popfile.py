import csv
import sys

sys.stderr = open(snakemake.log[0], "w")
snakemake.params.samplelist.loc[snakemake.params.inds].to_csv(snakemake.output.inds, sep="\t", quoting=csv.QUOTE_NONE, header = True, index = False)