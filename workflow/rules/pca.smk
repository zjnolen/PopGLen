rule pca_pcangsd:
	input:
		pruned_beagle=rules.merge_pruned_beagles.output.beagle
	output:
		cov=results + "/pcangsd/{population}.cov"
	log:
		logs + "/pcangsd/{population}.log"
	envmodules:
		# avail on rackham, should be replaced with conda/singul one day
		"bioinfo-tools",
		"PCAngsd/0.982"
	params:
		prefix=results + "/pcangsd/{population}"
	shell:
		"""
		pcangsd.py -b {input.pruned_beagle} -o {params.prefix} &> {log}
		"""