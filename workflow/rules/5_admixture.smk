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
		angsd_container
	params:
		prefix=results+"/analyses/ngsadmix/"+dataset+
			"_{population}{dp}_K{kvalue}",
		extra=config["params"]["ngsadmix"]["extra"],
		reps=config["params"]["ngsadmix"]["extra"]
	threads: 4
	resources:
		time=lambda wildcards, attempt: attempt*10080
	shell:
		"""
		export TMPDIR={resources.tmpdir}
		export reps={params.reps}
		workflow/scripts/ngsadmix.sh -likes {input.beagle} {params.extra} \
			-K {wildcards.kvalue} -P {threads} -o {params.prefix} 2> {log}
		"""

localrules: plot_admix

rule plot_admix:
	input:
		rules.ngsAdmix.output.qopt,
		rules.popfile.output.inds
	output:
		report(
			results+"/plots/ngsadmix/"+dataset+
				"_{population}{dp}_K{kvalue}.pdf",
			category="Admixture")
	conda:
		"../envs/r.yaml"
	script:
		"../scripts/plot_admix.R"

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