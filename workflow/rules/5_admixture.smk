rule ngsAdmix:
	input:
		beagle=rules.merge_pruned_beagles.output.beagle
	output:
		qopt=results+"/analyses/ngsadmix/"+dataset+
			"_{population}{dp}_K{kvalue}.qopt",
		fopt=results+"/analyses/ngsadmix/"+dataset+
			"_{population}{dp}_K{kvalue}.fopt.gz"
	log:
		logs+"/ngsadmix/"+dataset+"_{population}{dp}_K{kvalue}.log"
	container:
		"docker://zjnolen/angsd:0.937"
	params:
		prefix=results+"/analyses/ngsadmix/"+dataset+
			"_{population}{dp}_K{kvalue}"
	threads: 4
	resources:
		time=lambda wildcards, attempt: attempt*1440
	shell:
		"""
		export TMPDIR={resources.tmpdir}
		workflow/scripts/ngsadmix.sh -likes {input.beagle} \
			-K {wildcards.kvalue} -P {threads} -o {params.prefix} 2> {log}
		"""

rule evalAdmix:
	input:
		beagle=rules.merge_pruned_beagles.output.beagle,
		qopt=rules.ngsAdmix.output.qopt,
		fopt=rules.ngsAdmix.output.fopt
	output:
		results+"/analyses/ngsadmix/"+dataset+
			"_{population}{dp}_K{kvalue}.corres"
	log:
		logs+"/evaladmix/"+dataset+"_{population}{dp}_K{kvalue}.log"
	shell:
		"""
		~/polyommatini/working/shared/bin/evalAdmix -beagle {input.beagle} \
			-fname {input.fopt} -qname {input.qopt} -o {output} -P {threads} \
			2> {log}
		"""