rule ngsAdmix:
	input:
		beagle="results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_all{dp}_pruned.beagle.gz"
	output:
		qopt="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_all{dp}_K{kvalue}.qopt",
		fopt="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_all{dp}_K{kvalue}.fopt.gz"
	log:
		"logs/{dataset}/ngsadmix/{dataset}.{ref}_all{dp}_K{kvalue}.log"
	container:
		angsd_container
	params:
		prefix=lambda w, output: os.path.splitext(output.qopt)[0],
		extra=config["params"]["ngsadmix"]["extra"],
		reps=config["params"]["ngsadmix"]["reps"],
		minreps=config["params"]["ngsadmix"]["minreps"],
		thresh=config["params"]["ngsadmix"]["thresh"],
		conv=config["params"]["ngsadmix"]["conv"]
	threads: 4
	resources:
		time=lambda wildcards, attempt: attempt*2880
	shell:
		"""
		(export TMPDIR={resources.tmpdir}
		export reps={params.reps}
		export minreps={params.reps}
		export thresh={params.thresh}
		export conv={params.conv}
		workflow/scripts/ngsadmix.sh -likes {input.beagle} {params.extra} \
			-K {wildcards.kvalue} -P {threads} -o {params.prefix}) &> {log}
		"""

rule plot_admix:
	input:
		"results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_all{dp}_K{kvalue}.qopt",
		"results/datasets/{dataset}/poplists/{dataset}_all.indiv.list"
	output:
		report(
			"results/datasets/{dataset}/plots/ngsadmix/{dataset}.{ref}_all{dp}_K{kvalue}.svg",
			category="Admixture",
			subcategory="NGSadmix",
			labels={
				"Topic":"Admixture",
				"K-value":"{kvalue}",
				"Subsampling":"{dp}",
				"Type":"admix plot"
			})
	log:
		"logs/{dataset}/ngsadmix/{dataset}.{ref}_all{dp}_K{kvalue}_plot.log"
	conda:
		"../envs/r.yaml"
	script:
		"../scripts/plot_admix.R"

rule evalAdmix:
	input:
		beagle="results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_all{dp}_pruned.beagle.gz",
		qopt="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_all{dp}_K{kvalue}.qopt",
		fopt="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_all{dp}_K{kvalue}.fopt.gz"
	output:
		"results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_all{dp}_K{kvalue}.corres"
	log:
		"logs/{dataset}/evaladmix/{dataset}.{ref}_all{dp}_K{kvalue}.log"
	container:
		evaladmix_container
	shell:
		"""
		evalAdmix -beagle {input.beagle} -fname {input.fopt} \
			-qname {input.qopt} -o {output} -P {threads} &> {log}
		"""

rule plot_evalAdmix:
	input:
		corres="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_all{dp}_K{kvalue}.corres",
		qopt="results/datasets/{dataset}/analyses/ngsadmix/{dataset}.{ref}_all{dp}_K{kvalue}.qopt",
		pops="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list"
	output:
		report(
			"results/datasets/{dataset}/plots/evaladmix/{dataset}.{ref}_all{dp}_K{kvalue}_evaladmix.html",
			category="Admixture",
			subcategory="evalAdmix",
			labels={
				"Topic":"evalAdmix",
				"K-value":"{kvalue}",
				"Subsampling":"{dp}",
				"Type":"correlation matrix"
			})
	log:
		"logs/{dataset}/evaladmix/{dataset}.{ref}_all{dp}_K{kvalue}_plot.log"
	container:
		evaladmix_container
	script:
		"../scripts/plot_evaladmix.Rmd"