#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import sys
import math

parser = argparse.ArgumentParser(description='Calculates mean population sequencing depth from ANGSD .depthGlobal output file')

parser.add_argument("--input", help="The .depthGlobal output file from ANGSD to calculate the mean depth from.")
parser.add_argument("--min_mult", help="A multiplier to determine what proportion of the average to consider for minimum coverage.")
parser.add_argument("--max_mult", help="A multiplier to determine what multiple of the average to consider for maximum coverage.")
args = parser.parse_args()

df = pd.read_csv(sys.stdin, header = None, sep = "\t")
df = pd.DataFrame(df.sum())
df.columns = ['sites']
df["depth"] = df.index
df["scale_cov"] = df["sites"] * df["depth"]
avg_cov = np.nansum(df["scale_cov"]) / np.nansum(df["sites"])
min_cov = int(float(args.min_mult) * avg_cov)
max_cov = int(float(args.max_mult) * avg_cov)
avg_cov = float(avg_cov)

print(avg_cov, min_cov, max_cov, sep = "\t", file = sys.stdout)