
inds = get_samples_from_pop(wildcards.population)
samples.loc[inds].to_csv(output.inds, sep="\t", quoting=csv.QUOTE_NONE,
	header = True, index = False)