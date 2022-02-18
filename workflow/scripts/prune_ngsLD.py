#!/usr/bin/env python3
#
# Maintained by Zachary J. Nolen
#
# Copyright Zachary J. Nolen
#
# This script aims to efficiently prune SNPs from pairwise linkage
# disequilibrium measurements output by ngsLD. Runtime scales well with 
# dataset size, but is not multithreaded and memory usage seems to 
# require RAM equivalent to ~60-90% of uncompressed input file size


import datetime
import csv
import argparse
import pandas as pd
from graph_tool.all import *

parser = argparse.ArgumentParser(description='Prunes SNPs from ngsLD output to produce a list of sites in linkage equilibrium.')

parser.add_argument("--input", help="The .ld output file from ngsLD to be pruned.")
parser.add_argument("--output", help="The file to output pruned SNPs to.")
args = parser.parse_args()

begin_time = datetime.datetime.now()

# print("Reading in data...", flush=True)
# df = pd.read_table(args.input, header = None)
# df.columns = ['site1','site2','dist','r2pear','D','Dp','r2']
# df = df.drop(columns=['r2pear','D','Dp'])
# df['drop'] = (df['dist'] > 50000) | (df['r2'] < 0.1)

print("Reading in data...", flush=True)
G = load_graph_from_csv(args.input, directed = False, 
		eprop_types = ["int32_t","bool","bool","bool","double"],
		eprop_names = ["dist","na","na","na","r2"], hashed = True, 
		csv_options = {'delimiter': '\t'})

del G.ep["na"]

map_property_values(G.ep["r2"], G.ep["r2"], lambda x: abs(x))

drop_dist = G.new_edge_property("bool")
drop_r2 = G.new_edge_property("bool")
weight = G.new_vertex_property("double")

print("Filtering edges from graph...", flush=True)
map_property_values(G.ep["dist"], drop_dist, lambda x: x > 50000)
G.set_edge_filter(drop_dist, inverted=True)
G.purge_edges()
map_property_values(G.ep["r2"], drop_r2, lambda x: x < 0.1)
G.set_edge_filter(drop_r2, inverted=True)
G.purge_edges()
G.clear_filters()

# print("Making graph from dataset...", flush=True)
# name = G.add_edge_list(df.values, hashed=True, eprops=[dist,r2,drop])
# print("Filtering edges from graph...")
# G.set_edge_filter(drop,inverted=True)

print("Dropping neighbors from heaviest vertices...", flush=True)
while True:
	edges = G.num_edges()
	if edges == 0:
		break
	print(edges, flush = True)
	incident_edges_op(G, "out", "sum", G.ep["r2"], weight)
	max_weight = max(weight)
	print(max_weight, flush = True)
	heavy = find_vertex(G, weight, max_weight)
	heavy_neighbors = G.get_out_neighbors(heavy[0])
	G.remove_vertex(heavy_neighbors, fast = True)

print("Exporting kept sites to file...", flush=True)
pruned_df = pd.DataFrame([G.vp["name"][v] for v in G.get_vertices()])
pruned_df = pruned_df[0].str.split(pat=":", expand = True)
pruned_df.columns = ['chr','pos']
pruned_df.chr = pruned_df.chr.astype('string')
pruned_df.pos = pruned_df.pos.astype('int')
pruned_df = pruned_df.sort_values('pos')
pruned_df.to_csv(args.output, sep="\t", quoting=csv.QUOTE_NONE, 
	header = False, index = False)

print(datetime.datetime.now() - begin_time, flush=True)