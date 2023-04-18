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
		reps=config["params"]["ngsadmix"]["reps"]
	threads: 4
	resources:
		time=lambda wildcards, attempt: attempt*2880
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
				"_{population}{dp}_K{kvalue}.svg",
			category="Admixture",
			subcategory="NGSadmix",
			labels={
				"Topic":"Admixture",
				"K-value":"{kvalue}",
				"Subsampling":"{dp}",
				"Type":"admix plot"
			})
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
	container:
		evaladmix_container
	shell:
		"""
		evalAdmix -beagle {input.beagle} -fname {input.fopt} \
			-qname {input.qopt} -o {output} -P {threads} 2> {log}
		"""

localrules: plot_evalAdmix

rule plot_evalAdmix:
	input:
		corres=rules.evalAdmix.output,
		qopt=rules.ngsAdmix.output.qopt,
		pops=rules.popfile.output.inds
	output:
		report(
			results+"/plots/evaladmix/"+dataset+
				"_{population}{dp}_K{kvalue}_evaladmix.html",
			category="Admixture",
			subcategory="evalAdmix",
			labels={
				"Topic":"evalAdmix",
				"K-value":"{kvalue}",
				"Subsampling":"{dp}",
				"Type":"correlation matrix"
			})
	log:
		logs+"/evaladmix/"+dataset+"_{population}{dp}_K{kvalue}_plot.log"
	container:
		evaladmix_container
	script:
		"../scripts/plot_evaladmix.Rmd"