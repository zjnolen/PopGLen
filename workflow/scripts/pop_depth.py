#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import sys
import math

parser = argparse.ArgumentParser(description='Calculates mean population sequencing depth from ANGSD .depthGlobal output file')

parser.add_argument("--input", help="The .depthGlobal output file from ANGSD to calculate the mean depth from.")
args = parser.parse_args()

df = pd.read_csv(sys.stdin, header = None, sep = "\t")
df = pd.DataFrame(df.sum())
df.columns = ['sites']
df["depth"] = df.index
df["scale_cov"] = df["sites"] * df["depth"]
avg_cov = np.nansum(df["scale_cov"]) / np.nansum(df["sites"])
avg_cov = math.ceil(avg_cov)
print(avg_cov, file = sys.stdout)