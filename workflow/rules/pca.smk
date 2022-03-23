rule pca_pcangsd:
	input:
		pruned_beagle=rules.merge_pruned_beagles.output.beagle
	output:
		cov=results + "/pcangsd/{population}.cov"
	log:
		logs + "/pcangsd/{population}.log"
	params:
		prefix=results + "/pcangsd/{population}"
	shell:
		"""
		module load bioinfo-tools
		module load PCAngsd/0.982

		pcangsd.py -b {input.pruned_beagle} -o {params.prefix} &> {log}
		"""