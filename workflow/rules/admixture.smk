rule ngsAdmix:
	input:
		beagle=rules.merge_pruned_beagles.output.beagle
	output:
		qopt=results + "/ngsadmix/{population}_K{kvalue}.qopt"
	log:
		logs + "/ngsadmix/{population}_K{kvalue}.log"
	conda:
		"../envs/angsd.yaml"
	params:
		prefix=results + "/ngsadmix/{population}_K{kvalue}"
	shell:
		"""
		workflow/scripts/wrapper_ngsAdmix.sh -likes {input.beagle} \
			-K ${wildcards.kvalue} -P {threads} -o {params.prefix}
		"""