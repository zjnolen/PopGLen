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
	threads: 4
	resources:
		time="24:00:00"
	shell:
		"""
		workflow/scripts/wrapper_ngsAdmix.sh -likes {input.beagle} \
			-K {wildcards.kvalue} -P {threads} -o {params.prefix}
		"""